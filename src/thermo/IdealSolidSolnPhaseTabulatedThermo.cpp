/**
 *  @file IdealSolidSolnPhaseTabulatedThermo.cpp
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

#include "cantera/thermo/IdealSolidSolnPhaseTabulatedThermo.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/vec_functions.h"
#include "cantera/base/ctml.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/GeneralSpeciesThermo.h"

namespace Cantera
{
IdealSolidSolnPhaseTabulatedThermo::IdealSolidSolnPhaseTabulatedThermo(const IdealSolidSolnPhaseTabulatedThermo& right)
{
	*this = right;
}

IdealSolidSolnPhaseTabulatedThermo& IdealSolidSolnPhaseTabulatedThermo::operator=(const IdealSolidSolnPhaseTabulatedThermo& right)
{
	if (&right == this) {
		return *this;
	}

	m_h0_RT         = right.m_h0_RT;
	m_cp0_R         = right.m_cp0_R;
	m_g0_RT         = right.m_g0_RT;
	m_s0_R          = right.m_s0_R;
	m_pp            = right.m_pp;
	m_xlast			= right.m_xlast;

	return *this;

}

ThermoPhase* IdealSolidSolnPhaseTabulatedThermo::duplMyselfAsThermoPhase() const
{
	return new IdealSolidSolnPhaseTabulatedThermo(*this);
}

void IdealSolidSolnPhaseTabulatedThermo::getPartialMolarEnthalpies(doublereal* result) const {
	doublereal vdp = (pressure() - m_spthermo->refPressure()) / molarDensity();
	doublereal rt = temperature() * GasConstant;
	const vector_fp& h_RT = enthalpy_RT_ref();
	for (size_t k = 0; k < m_kk; k++) {
		result[k] = rt*h_RT[k] + vdp;
	}
}

void IdealSolidSolnPhaseTabulatedThermo::getPartialMolarEntropies(doublereal* result) const {
	doublereal xx;
	doublereal r = GasConstant;
	const vector_fp& s_R = entropy_R_ref();
	for (size_t k = 0; k < m_kk; k++) {
		xx = moleFraction(k);
		if(xx < SmallNumber) xx = SmallNumber;
		result[k] = r*(s_R[k] - log(xx));
	}
}

void IdealSolidSolnPhaseTabulatedThermo::setMoleFractions(const doublereal* const x)
{
	Phase::setMoleFractions(x);
	_updateThermo();
	calcDensity();
}

void IdealSolidSolnPhaseTabulatedThermo::setMoleFractions_NoNorm(const doublereal* const x)
{
	Phase::setMoleFractions_NoNorm(x);
	_updateThermo();
	calcDensity();
}

void IdealSolidSolnPhaseTabulatedThermo::setMassFractions(const doublereal* const y)
{
	Phase::setMassFractions(y);
	_updateThermo();
	calcDensity();
}

void IdealSolidSolnPhaseTabulatedThermo::setMassFractions_NoNorm(const doublereal* const y)
{
	Phase::setMassFractions_NoNorm(y);
	_updateThermo();
	calcDensity();
}

void IdealSolidSolnPhaseTabulatedThermo::setConcentrations(const doublereal* const conc)
{
	Phase::setConcentrations(conc);
	_updateThermo();
	calcDensity();
}

void IdealSolidSolnPhaseTabulatedThermo::_updateThermo() const
{
	doublereal tnow = temperature();
	doublereal xnow = moleFraction(m_kk_mod);
	doublereal c[4];
	doublereal dS_corr = 0.0;
	if (m_tlast != tnow || m_xlast != xnow) {
		GeneralSpeciesThermo* genSpthermo = dynamic_cast<GeneralSpeciesThermo*>(m_spthermo);
		if (!genSpthermo) {
			throw CanteraError("IdealSolidSolnPhaseTabulatedThermo::_updateThermo",
					"unable to modify IdealSolidSolnPhaseTabulatedThermo thermo");
		} else {
			c[0]=tnow;
			c[1]=interp_h(xnow) * 1e3; // 1e3 for conversion J/mol -> J/kmol
			if (xnow==0) {
				dS_corr = -BigNumber;
			} else if (xnow==1) {
				dS_corr = BigNumber;
			} else {
				dS_corr = GasConstant*std::log(xnow/(1.0-xnow));
			}
			c[2]=interp_s(xnow) * 1e3 + dS_corr; // 1e3 for conversion J/K/mol -> J/K/kmol
			c[3]=0.0;
			genSpthermo->modifyParams(m_kk_mod, &c[0]);
			m_spthermo->update(tnow, &m_cp0_R[0], &m_h0_RT[0],
					&m_s0_R[0]);
			for (size_t k = 0; k < m_kk; k++) {
				m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
			}
			m_xlast = xnow;
			m_tlast = tnow;
		}
	}
}

void IdealSolidSolnPhaseTabulatedThermo::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
	size_t offset;
	doublereal a=0,b=0,c=0;
	std::ifstream t;
	std::string line;
	int count = 0;

	if (id_.size() > 0) {
		if (phaseNode.id() != id_) {
			throw CanteraError("IdealSolidSolnPhaseTabulatedThermo::initThermoXML",
					"phasenode and Id are incompatible");
		}
	}
	if (phaseNode.hasChild("thermo")) {
		XML_Node& thermoNode = phaseNode.child("thermo");
		std::string mString = thermoNode["model"];
		if (lowercase(mString) != "idealsolidsolutiontabulatedthermo") {
			throw CanteraError("IdealSolidSolnPhaseTabulatedThermo::initThermoXML",
					"Unknown thermo model: " + mString);
		}
		if (thermoNode.hasChild("data")) {
			XML_Node& scNode = thermoNode.child("data");
			std::vector<std::string> nameFile;
			getStringArray(scNode, nameFile);
			if (nameFile.size() != 1) {
				throw CanteraError("IdealSolidSolnPhaseTabulatedThermo::initThermoXML",
						"badly formed data XML node");
			}
			m_dataFile = nameFile[0];
		} else {
			throw CanteraError("initThermoXML",
					"Unspecified data file");
		}
		if (thermoNode.hasChild("modifiable_species")) {
			std::string modifiable_species_name = thermoNode.child("modifiable_species").value();
			m_kk_mod = speciesIndex(modifiable_species_name);
			if (m_kk_mod == npos) {
				throw CanteraError("IdealSolidSolnPhaseTabulatedThermo::initThermoXML",
						"Species " + modifiable_species_name + " not found.");
			}
		} else {
			throw CanteraError("IdealSolidSolnPhaseTabulatedThermo::initThermoXML",
					"Unspecified modifiable species");
		}
	} else {
		throw CanteraError("IdealSolidSolnPhaseTabulatedThermo::initThermoXML",
				"Unspecified thermo model");
	}

	/*
	 * Form of the standard concentrations. Must have one of:
	 *
	 *     <standardConc model="unity" />
	 *     <standardConc model="molar_volume" />
	 *     <standardConc model="solvent_volume" />
	 */
	if (phaseNode.hasChild("standardConc")) {
		XML_Node& scNode = phaseNode.child("standardConc");
		std::string formStringa = scNode.attrib("model");
		std::string formString = lowercase(formStringa);
		if (formString == "unity") {
			m_formGC = 0;
		} else if (formString == "molar_volume") {
			m_formGC = 1;
		} else if (formString == "solvent_volume") {
			m_formGC = 2;
		} else {
			throw CanteraError("IdealSolidSolnPhase::initThermoXML",
					"Unknown standardConc model: " + formStringa);
		}
	} else {
		throw CanteraError("IdealSolidSolnPhase::initThermoXML",
				"Unspecified standardConc model");
	}

	/*
	 * Initialize all of the lengths now that we know how many species
	 * there are in the phase.
	 */
	initLengths();

	/*
	 * Now go get the molar volumes
	 */
	XML_Node& speciesList = phaseNode.child("speciesArray");
	XML_Node* speciesDB = get_XML_NameID("speciesData", speciesList["datasrc"],
			&phaseNode.root());

	for (size_t k = 0; k < m_kk; k++) {
		XML_Node* s =  speciesDB->findByAttr("name", speciesName(k));
		XML_Node* ss = s->findByName("standardState");
		m_speciesMolarVolume[k] = getFloat(*ss, "molarVolume", "toSI");
	}

	/*
	 * Read data file and create the interpolation table
	 */
	t.open(m_dataFile.c_str());
	if(t.fail())
	{
		throw CanteraError("IdealSolidSolnPhase::initThermoXML",
				"Unknown data file: " + m_dataFile);
	}
	while (getline(t,line))
	{
		if((offset = line.find("*", 0)) != std::string::npos)
		{
			continue;
		}
		else
		{
			count = std::sscanf(line.c_str(), "%lf %lf %lf", &a, &b, &c);
			if (count==3)
			{
				molefrac_h.push_back(std::make_pair(a,b));
				molefrac_s.push_back(std::make_pair(a,c));
			}
			else
			{
				throw CanteraError("IdealSolidSolnPhase::initThermoXML",
						"Data format error: " + m_dataFile);
			}
		}
	}
	// Sort the xLi, h, s data in the order of increasing xLi
	std::sort(molefrac_h.begin(), molefrac_h.end());
	std::sort(molefrac_s.begin(), molefrac_s.end());
	t.close();

    /*
     * Call the base initThermo, which handles setting the initial
     * state.
     */
	ThermoPhase::initThermoXML(phaseNode, id_);
}

doublereal IdealSolidSolnPhaseTabulatedThermo::interp_h(doublereal x) const
{
	// Assumes that "table" is sorted by .first
	// Check if x is out of bound
	if (x > molefrac_h.back().first) return molefrac_h.back().second;
	if (x < molefrac_h[0].first) return molefrac_h[0].second;
	std::vector<std::pair<double, double> >::const_iterator it, it2;
	// INFINITY is defined in math.h in the glibc implementation
	it = std::lower_bound(molefrac_h.begin(), molefrac_h.end(), std::make_pair(x, molefrac_h[0].second));
	// Corner case
	if (it == molefrac_h.begin()) return it->second;
	it2 = it;
	--it2;
	return it2->second + (it->second - it2->second)*(x - it2->first)/(it->first - it2->first);
}

doublereal IdealSolidSolnPhaseTabulatedThermo::interp_s(doublereal x) const
{
	// Assumes that "table" is sorted by .first
	// Check if x is out of bound
	if (x > molefrac_s.back().first) return molefrac_s.back().second;
	if (x < molefrac_s[0].first) return molefrac_s[0].second;
	std::vector<std::pair<double, double> >::const_iterator it, it2;
	// INFINITY is defined in math.h in the glibc implementation
	it = std::lower_bound(molefrac_s.begin(), molefrac_s.end(), std::make_pair(x, molefrac_s[0].second));
	// Corner case
	if (it == molefrac_s.begin()) return it->second;
	it2 = it;
	--it2;
	return it2->second + (it->second - it2->second)*(x - it2->first)/(it->first - it2->first);
}

void IdealSolidSolnPhaseTabulatedThermo::initLengths()
{
	/*
	 * Obtain the reference pressure by calling the ThermoPhase
	 * function refPressure, which in turn calls the
	 * species thermo reference pressure function of the
	 * same name.
	 */
	m_Pref = refPressure();
	m_h0_RT.resize(m_kk);
	m_g0_RT.resize(m_kk);
	m_expg0_RT.resize(m_kk);
	m_cp0_R.resize(m_kk);
	m_s0_R.resize(m_kk);
	m_pe.resize(m_kk, 0.0);
	m_pp.resize(m_kk);
	m_speciesMolarVolume.resize(m_kk);
}

}
