/**
 *  @file misc.cpp
 */

/*
 * $Author$
 * $Revision$
 * $Date$
 */

#ifdef WIN32
#pragma warning(disable:4786)
#pragma warning(disable:4503)
#endif

#include "global.h"
#include "ctexceptions.h"
#include "stringUtils.h"
#include "units.h"
#include "xml.h"
#include "ctml.h"
#include "SpeciesThermoFactory.h"
#include "ThermoFactory.h"
#include "FalloffFactory.h"

//#ifndef WIN32
//#include "ctdir.h"
//#endif

#include <fstream>
using namespace std;

namespace Cantera {

    /**
     * Class to hold global data. Class Application is the top-level
     * class that stores data that should persist for the duration of
     * the process. The class should not be instantiated directly;
     * instead, it is instantiated as needed by the functions declared
     * here. At most one instance is created, and it is not destroyed
     * until the process terminates.
     */
    class Application {
    public:
        Application() : linelen(0), stop_on_error(false),
#ifdef WIN32
                        tmp_dir(".") 
#else
            tmp_dir("/tmp") 
#endif
            {
                char* tmpdir = getenv("TMP");
                if (tmpdir == 0) 
                    tmpdir = getenv("TEMP");
                if (tmpdir != 0)
                    tmp_dir = string(tmpdir);
            }

        virtual ~Application() {
            map<string, XML_Node*>::iterator pos;
            for (pos = xmlfiles.begin(); pos != xmlfiles.end(); ++pos) {
                pos->second->unlock();
                delete pos->second;
                pos->second = 0;
            }
        }
        vector<string> inputDirs;
        vector<string> errorMessage;
        vector<string> warning;
        vector<string> errorRoutine;
        string msglog;
        size_t linelen;
        bool stop_on_error;
        map<string, string>     options;
        string tmp_dir;
        map<string, XML_Node*> xmlfiles;
    };
        
            
    /// Returns a pointer to the one and only instance of Application
    Application* app();

    void setDefaultDirectories();

    static Application* __app = 0;
    Unit* Unit::__u = 0;

    static void appinit() {
        if (__app == 0) __app = new Application;
    }

    /**
     * This function deletes the global information. It should be called
     * at the end of the application, especially if leak checking is
     * to be done.
     */
    void appdelete() {
        if (__app) {
            delete __app;
            __app = 0;
        }
        SpeciesThermoFactory::deleteFactory();
        ThermoFactory::deleteFactory();
        FalloffFactory::deleteFalloffFactory();
        Unit::deleteUnit();
    }

    Application* app() {
        if (__app == 0) {
            __app = new Application;
            setDefaultDirectories();
        }
        return __app;
    }

    XML_Node* get_XML_File(string file) {
        string path = findInputFile(file);
        string ff = path;
        if (app()->xmlfiles.find(path) 
            == app()->xmlfiles.end()) {
            /*
             * Check whether or not the file is XML. If not, it will
             * be first processed with the preprocessor. We determine
             * whether it is an XML file by looking at the file extension.
             */
            string::size_type idot = path.rfind('.');
            string ext;
            if (idot != string::npos) {
                ext = path.substr(idot, path.size());
            } else {
                ext = "";
                idot = path.size();
            }
            if (ext != ".xml" && ext != ".ctml") {
                /*
                 * We will assume that we are trying to open a cti file.
                 * First, determine the name of the xml file, ff, derived from
                 * the cti file.
                 * In all cases, we will write the xml file to the current
                 * directory.
                 */
                string::size_type islash = path.rfind('/');
                if (islash != string::npos) 
                    ff = string("./")+path.substr(islash+1,idot-islash - 1) + ".xml";
                else {
                    ff = string("./")+path.substr(0,idot) + ".xml";
                }
#ifdef DEBUG_PATHS
                cout << "get_XML_File(): Expected location of xml file = "
                     << ff << endl;
#endif
                /*
                 * Do a search of the existing XML trees to determine if we have
                 * already processed this file. If we have, return a pointer to
                 * the processed xml tree.
                 */
                if (app()->xmlfiles.find(ff) != app()->xmlfiles.end()) {
#ifdef DEBUG_PATHS
                    cout << "get_XML_File(): File, " << ff << ", was previously read."
                         << " Retrieving the storred xml tree." << endl;
#endif
                    return __app->xmlfiles[ff];	  
                }
                /*
                 * Ok, we didn't find the processed XML tree. Do the conversion
                 * to xml, possibly overwriting the file, ff, in the process.
                 */
                ctml::ct2ctml(path.c_str());
            }
            else {
                ff = path;
            }
            /*
             * Take the XML file ff, open it, and process it, creating an
             * XML tree, and then adding an entry in the map. We will store
             * the absolute pathname as the key for this map.
             */
            ifstream s(ff.c_str());
            XML_Node* x = new XML_Node("doc");
            if (s) {
                x->build(s);
                x->lock();
                __app->xmlfiles[ff] = x;
            }
            else {
                string estring = "cannot open "+ff+" for reading.";
                estring += "Note, this error indicates a possible configuration problem."; 
                throw CanteraError("get_XML_File", estring);
            }
        }
        /*
         * Return the XML node pointer. At this point, we are sure that the
         * lookup operation in the return statement will return a valid
         * pointer. 
         */
        return __app->xmlfiles[ff];
    }

    void close_XML_File(string file) {
        if (file == "all") {
            map<string, XML_Node*>::iterator 
                b = app()->xmlfiles.begin(), e = app()->xmlfiles.end();
            for(; b != e; ++b) {
                b->second->unlock();
                delete b->second;
                __app->xmlfiles.erase(b->first);
            }
        }
        else if (app()->xmlfiles.find(file) 
            != app()->xmlfiles.end()) {
            __app->xmlfiles[file]->unlock();
            delete __app->xmlfiles[file];
            __app->xmlfiles.erase(file);
        }
    }

    void setTmpDir(string tmp) { app()->tmp_dir = tmp; }
    string tmpDir() { appinit(); return app()->tmp_dir; }

    int nErrors() {
        return static_cast<int>(app()->errorMessage.size());
    }

    void popError() {
        appinit();
        if (nErrors() > 0) {
            __app->errorMessage.pop_back();
            __app->errorRoutine.pop_back();
        }
    }

    string lastErrorMessage() {
        appinit();
        if (nErrors() > 0) {
            string head = 
                "\n\n************************************************\n"
                "                Cantera Error!                  \n"
                "************************************************\n\n";
            return head+string("\nProcedure: ")+__app->errorRoutine.back()
                +string("\nError:   ")+__app->errorMessage.back();
        }
        else
            return "<no Cantera error>";
    }

    void showErrors(ostream& f) {
        appinit(); 
        int i = static_cast<int>(__app->errorMessage.size());
        if (i == 0) return;
        f << endl << endl;
        f << "************************************************" << endl;
        f << "                Cantera Error!                  " << endl;
        f << "************************************************" << endl << endl;
        int j;
        for (j = 0; j < i; j++) {
            f << endl;
            f << "Procedure: " << __app->errorRoutine[j] << endl;
            f << "Error:     " << __app->errorMessage[j] << endl;
        } 
        f << endl << endl;
        __app->errorMessage.clear();
        __app->errorRoutine.clear();
    }

    void showErrors() {
        appinit(); 
        int i = static_cast<int>(__app->errorMessage.size());
        if (i == 0) return;
        writelog("\n\n");
        writelog("************************************************\n");
        writelog("                Cantera Error!                  \n");
        writelog("************************************************\n\n");
        int j;
        for (j = 0; j < i; j++) {
            writelog("\n");
            writelog(string("Procedure: ")+ __app->errorRoutine[j]+" \n");
            writelog(string("Error:     ")+__app->errorMessage[j]+" \n");
        } 
        writelog("\n\n");
        __app->errorMessage.clear();
        __app->errorRoutine.clear();
    }

    void setError(string r, string msg) {
        appinit();
        __app->errorMessage.push_back(msg);
        __app->errorRoutine.push_back(r);
    }


    /**
     * Set the default directories for input data files. 
     */
    void setDefaultDirectories() {
        appinit();
        vector<string>& dirs = __app->inputDirs;

        // always look in the local directory first
        dirs.push_back(".");


#ifdef WIN32
        //
        // Under Windows, the Cantera setup utility puts data files in
        // a directory 'Cantera\data' below the one the environment
        // variable COMMONPROGRAMFILES points to. (This is usually
        // C:\Program Files\Common Files.) If this environment
        // variable is defined, then this directory is assumed to
        // exist and is added to the search path.
        //
        const char* comfiles = getenv("COMMONPROGRAMFILES");
        if (comfiles != 0) {
            string cfiles = string(comfiles);

            // remove quotes if necessary
            if (cfiles[0] == '\'') 
                cfiles = cfiles.substr(1,1000);
            if (cfiles[cfiles.size()-1] == '\'') cfiles[cfiles.size()-1] = '\n';

            string datadir = string(comfiles) + "/Cantera/data";
            string tmpldir = string(comfiles) + "/Cantera/templates";
            dirs.push_back(datadir);
            dirs.push_back(tmpldir);
        }
#endif

#ifdef DARWIN
        //
        // add a default data location for Mac OS X
        //
        if (DARWIN == 1)
            dirs.push_back("/Applications/Cantera/data");
#endif

        //
        // if environment variable CANTERA_DATA is defined, then add
        // it to the search path
        //
        if (getenv("CANTERA_DATA") != 0) {
            string datadir = string(getenv("CANTERA_DATA"));
            dirs.push_back(datadir);
        }

        // CANTERA_DATA is defined in file config.h. This file is written
        // during the build process (unix), and points to the directory
        // specified by the 'prefix' option to 'configure', or else to
        // /usr/local/cantera. 
#ifdef CANTERA_DATA
        string datadir = string(CANTERA_DATA);
        dirs.push_back(datadir);
#endif
    }




    void addDirectory(string dir) {
        appinit();
        if (__app->inputDirs.size() == 0) setDefaultDirectories();
        __app->inputDirs.push_back(stripnonprint(dir));
    }

    /**
     *  findInputFile():
     *
     *    This routine will search for a file in the default
     *    locations specified for the application.
     *    See the routine setDefaultDirectories() listed above.
     *
     *    The default set of directories specified for the application
     *    will be searched if a '/' or an '\\' is not found in
     *    name. If either is found then a relative path name is
     *    presumed and the default directories are not searched.
     *
     *    The presence of the file is determined by whether the file
     *    can be opened for reading by the current user.
     *
     *    Return
     *    -------
     *      The absolute path name of the first matching
     *      file is returned. If a relative path name
     *      is indicated, the relative path name is returned.
     *  
     *      If the file is not found, a message is written to 
     *      stdout and  a CanteraError exception is thrown.
     */
    string findInputFile(string name) {
        appinit();
        string::size_type islash = name.find('/');
        string::size_type ibslash = name.find('\\');
        string inname;
        vector<string>& dirs = __app->inputDirs;
        if (dirs.size() == 0) setDefaultDirectories();

        int nd;
        if (islash == string::npos && ibslash == string::npos) {
            nd = static_cast<int>(dirs.size());
            int i;
            inname = "";
            for (i = 0; i < nd; i++) {
                inname = dirs[i] + "/" + name;
                ifstream fin(inname.c_str());
                if (fin) {
                    fin.close();
                    return inname;
                }
            }
            string msg;
            msg = "\nInput file " + name 
                  + " not found in director";
            msg += (nd == 1 ? "y " : "ies ");
            for (i = 0; i < nd; i++) {
                msg += "\n'" + dirs[i] + "'";
                if (i < nd-1) msg += ", ";
            }
            msg += "\n\n";
            msg += "To fix this problem, either:\n";
            msg += "    a) move the missing files into the local directory;\n";
            msg += "    b) define environment variable CANTERA_DATA to\n";
            msg += "         point to the directory containing the file.";
            throw CanteraError("findInputFile", msg);
            return "";
        }
        else {
            return name;
        }
    }

    void writelog(const char* msg) {writelog(string(msg));}

    doublereal toSI(string unit) {
        doublereal f = Unit::units()->toSI(unit);
        if (f) return f;
        else return 1.0;
    }

    doublereal actEnergyToSI(string unit) {
        doublereal f = Unit::units()->actEnergyToSI(unit);
        if (f) return f;
        else return 1.0;
    }


    string canteraRoot() {
        char* ctroot = 0;
        ctroot = getenv("CANTERA_ROOT");
        if (ctroot != 0) { return string(ctroot); }
        else {
#ifdef CANTERA_ROOT
            return string(CANTERA_ROOT);
#else
            return "";
#endif
        }
    }


    /// exceptions

    CanteraError::CanteraError(string proc, string msg) {
        setError(proc, msg);
    }
    
    ArraySizeError::ArraySizeError(string proc, int sz, int reqd) :
        CanteraError(proc, "Array size ("+int2str(sz)+
            ") too small. Must be at least "+int2str(reqd)) {}

    ElementRangeError::ElementRangeError(string func, int m, int mmax) :
        CanteraError(func, "Element index " + int2str(m) + 
            " outside valid range of 0 to " + int2str(mmax-1)) {}
}


