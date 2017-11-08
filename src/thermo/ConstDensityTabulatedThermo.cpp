/**
 *  @file ConstDensityTabulatedThermo.cpp
 * Definition file for a derived class of ThermoPhase that assumes either
 * an ideal gas or ideal solution approximation and handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties (see \ref thermoprops and
 * class \link Cantera::ConstDensityTabulatedThermo ConstDensityTabulatedThermo\endlink).
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

#include "cantera/thermo/ConstDensityTabulatedThermo.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/ctml.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/MultiSpeciesThermo.h"
#include <memory>

namespace Cantera
{

ConstDensityTabulatedThermo::ConstDensityTabulatedThermo(const ConstDensityTabulatedThermo& right)
{
	*this = right;
}

ConstDensityTabulatedThermo& ConstDensityTabulatedThermo::operator=(const ConstDensityTabulatedThermo& right)
{
	if (&right == this) {
		return *this;
	}

	m_h0_RT         = right.m_h0_RT;
	m_cp0_R         = right.m_cp0_R;
	m_g0_RT         = right.m_g0_RT;
	m_s0_R          = right.m_s0_R;
	m_xlast			= right.m_xlast;

	return *this;
}

ThermoPhase* ConstDensityTabulatedThermo::duplMyselfAsThermoPhase() const
{
	return new ConstDensityTabulatedThermo(*this);
}

void ConstDensityTabulatedThermo::getActivityConcentrations(doublereal* c) const
{
	getMoleFractions(c);
}

void ConstDensityTabulatedThermo::getActivityCoefficients(doublereal* ac) const
{
	for (size_t k = 0; k < m_kk; k++) {
		ac[k] = 1.0;
	}
}

void ConstDensityTabulatedThermo::getPartialMolarEnthalpies(doublereal* result) const {
	doublereal vdp = (pressure() - m_spthermo->refPressure()) / molarDensity();
	doublereal rt = temperature() * GasConstant;
	const vector_fp& h_RT = enthalpy_RT();
	for (size_t k = 0; k < m_kk; k++) {
		result[k] = rt*h_RT[k] + vdp;
	}
}

void ConstDensityTabulatedThermo::getPartialMolarEntropies(doublereal* result) const {
	doublereal xx;
	doublereal r = GasConstant;
	const vector_fp& s_R = entropy_R();
	for (size_t k = 0; k < m_kk; k++) {
		xx = moleFraction(k);
		if(xx < SmallNumber) xx = SmallNumber;
		result[k] = r*(s_R[k] - log(xx));
	}
}

doublereal ConstDensityTabulatedThermo::standardConcentration(size_t k) const
{
	return 1;
}

void ConstDensityTabulatedThermo::setMoleFractions(const doublereal* const x)
{
	Phase::setMoleFractions(x);
	_updateThermo();
}

void ConstDensityTabulatedThermo::setMoleFractions_NoNorm(const doublereal* const x)
{
	Phase::setMoleFractions_NoNorm(x);
	_updateThermo();
}

void ConstDensityTabulatedThermo::setMassFractions(const doublereal* const y)
{
	Phase::setMassFractions(y);
	_updateThermo();
}

void ConstDensityTabulatedThermo::setMassFractions_NoNorm(const doublereal* const y)
{
	Phase::setMassFractions_NoNorm(y);
	_updateThermo();
}

void ConstDensityTabulatedThermo::setConcentrations(const doublereal* const conc)
{
	Phase::setConcentrations(conc);
	_updateThermo();
}

void ConstDensityTabulatedThermo::_updateThermo() const
{
	doublereal tnow = temperature();
	doublereal xnow = moleFraction(m_kk_mod);
	doublereal c[4];
	doublereal dS_corr = 0.0;
	doublereal tlow = 0.0, thigh = 0.0;
	int type = 0;
	if (m_tlast != tnow || m_xlast != xnow) {
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
		type = m_spthermo->reportType(m_kk_mod);
		tlow = m_spthermo->minTemp(m_kk_mod);
		thigh = m_spthermo->maxTemp(m_kk_mod);
		shared_ptr<SpeciesThermoInterpType> stit(
				newSpeciesThermoInterpType(type, tlow, thigh, OneAtm, c));
		m_spthermo->modifySpecies(m_kk_mod, stit);
		m_spthermo->update(tnow, &m_cp0_R[0], &m_h0_RT[0],
				&m_s0_R[0]);
		for (size_t k = 0; k < m_kk; k++) {
			m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
		}
		m_xlast = xnow;
		m_tlast = tnow;
	}
}

void ConstDensityTabulatedThermo::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
	ConstDensityTabulatedThermo::initThermo();
	size_t offset;
	doublereal a=0,b=0,c=0;
	std::ifstream t;
	std::string line;
	int count = 0;

	if (phaseNode.hasChild("thermo")) {
		XML_Node& thermoNode = phaseNode.child("thermo");
		std::string model = thermoNode["model"];
		if (model == "ConstDensityTabulatedThermo") {
			if (thermoNode.hasChild("data")) {
				XML_Node& scNode = thermoNode.child("data");
				std::vector<std::string> nameFile;
				getStringArray(scNode, nameFile);
				if (nameFile.size() != 1) {
					throw CanteraError("ConstDensityTabulatedThermo::initThermoXML",
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
					throw CanteraError("ConstDensityTabulatedThermo::initThermoXML",
							"Species " + modifiable_species_name + " not found.");
				}
			} else {
				throw CanteraError("ConstDensityTabulatedThermo::initThermoXML",
						"Unspecified modifiable species");
			}
		} else {
			throw CanteraError("ConstDensityTabulatedThermo::initThermoXML",
					"Unknown thermo model : " + model);
		}
	}
	t.open(m_dataFile.c_str());
	if(t.fail())
	{
		throw CanteraError("initThermoXML",
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
				throw CanteraError("initThermoXML",
						"Data format error: " + m_dataFile);
			}
		}
	}
	// Sort the xLi, h, s data in the order of increasing xLi
	std::sort(molefrac_h.begin(), molefrac_h.end());
	std::sort(molefrac_s.begin(), molefrac_s.end());
	t.close();
	ThermoPhase::initThermoXML(phaseNode, id_);
}

doublereal ConstDensityTabulatedThermo::interp_h(doublereal x) const
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

doublereal ConstDensityTabulatedThermo::interp_s(doublereal x) const
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

void ConstDensityTabulatedThermo::setParametersFromXML(const XML_Node& eosdata)
{
	eosdata._require("model","ConstDensityTabulatedThermo");
	doublereal rho = getFloat(eosdata, "density", "toSI");
	setDensity(rho);
}

}
