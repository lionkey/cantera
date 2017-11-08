/*
 * IdealSolidSolnPhaseTabulatedThermo.h
 *
 *  Created on: 04.06.2016
 *      Author: mmayur
 */

/**
 *  @file IdealSolidSolnPhaseTabulatedThermo.h
 * Definition file for a derived class of ThermoPhase that assumes either
 * an ideal gas or ideal solution approximation and handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties (see \ref thermoprops and
 * class \link Cantera::IdealSolidSolnPhaseTabulatedThermo IdealSolidSolnPhaseTabulatedThermo\endlink).
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */
#ifndef CT_IDEALSOLIDSOLNPHASETABULATEDTHERMO_H
#define CT_IDEALSOLIDSOLNPHASETABULATEDTHERMO_H

#include "IdealSolidSolnPhase.h"
#include "cantera/base/utilities.h"
#include <vector>
#include <iterator>
#include <iostream>
#include <fstream>

namespace Cantera
{

//!   Overloads the virtual methods of class ThermoPhase to implement the
//!   incompressible equation of state.
/**
 * <b> Specification of Solution Thermodynamic Properties </b>
 *
 * The density is assumed to be constant, no matter what the concentration of the solution.
 *
 * @ingroup thermoprops
 */
class IdealSolidSolnPhaseTabulatedThermo : public IdealSolidSolnPhase
{
public:
    //! Constructor.
	IdealSolidSolnPhaseTabulatedThermo() {}

    //! Copy Constructor
    /*!
     * @param right Object to be copied
     */
	IdealSolidSolnPhaseTabulatedThermo(const IdealSolidSolnPhaseTabulatedThermo& right);

    //! Assignment Operator
    /*!
     * @param right Object to be copied
     */
	IdealSolidSolnPhaseTabulatedThermo& operator=(const IdealSolidSolnPhaseTabulatedThermo& right);

    //! Duplication routine for objects which inherit from ThermoPhase
    /*!
     *  This virtual routine can be used to duplicate objects
     *  derived from ThermoPhase even if the application only has
     *  a pointer to ThermoPhase to work with.
     */
    virtual ThermoPhase* duplMyselfAsThermoPhase() const;

    //! Was added because of the DENIS LFP Runaway simulation.
    virtual void getPartialMolarEnthalpies(doublereal* result) const;

    //! Was added because of the DENIS LFP Runaway simulation.
    virtual void getPartialMolarEntropies(doublereal* result) const;

    //! Initialize the ThermoPhase object after all species have been set up
    /*!
     * @internal Initialize.
     *
     * This method is provided to allow
     * subclasses to perform any initialization required after all
     * species have been added. For example, it might be used to
     * resize internal work arrays that must have an entry for
     * each species.  The base class implementation does nothing,
     * and subclasses that do not require initialization do not
     * need to overload this method.  When importing a CTML phase
     * description, this method is called from ThermoPhase::initThermoXML(),
     * which is called from importPhase(),
     * just prior to returning from function importPhase().
     */

    virtual void initThermoXML(XML_Node& phaseNode, const std::string& id_);
    virtual void setMoleFractions(const doublereal* const x);
    virtual void setMoleFractions_NoNorm(const doublereal* const x);
    virtual void setMassFractions(const doublereal* const y);
    virtual void setMassFractions_NoNorm(const doublereal* const y);
    virtual void setConcentrations(const doublereal* const conc);

    //! Species thermodynamics interpolation functions
	doublereal interp_h(doublereal x) const;
	doublereal interp_s(doublereal x) const;

    //! Current modifiable species index
    size_t m_kk_mod;

protected:

    //! Current modifiable species mole fraction
    mutable doublereal m_xlast;

    //! File name for IdealSolidSolnPhaseTabulatedThermo thermo data
    std::string m_dataFile;

    //! Vector pairs for IdealSolidSolnPhaseTabulatedThermo thermo
    std::vector<std::pair<doublereal,doublereal> > molefrac_h;
    std::vector<std::pair<doublereal,doublereal> > molefrac_s;

    //! Function to update the reference state thermo functions
    //! Changed from private to protected so that the methods of the parent class could use this overridden function
    void _updateThermo() const;

private:

    //! This internal function adjusts the lengths of arrays
    void initLengths();

};
}

//#endif

#endif /* INCLUDE_CANTERA_THERMO_IDEALSOLIDSOLNPHASETABULATEDTHERMO_H_ */
