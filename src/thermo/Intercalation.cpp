/**
 *  @file Intercalation.cpp
 * Definition file for a derived class of ThermoPhase that assumes either
 * an ideal gas or ideal solution approximation and handles
 * variable pressure standard state methods for calculating
 * thermodynamic properties (see \ref thermoprops and
 * class \link Cantera::Intercalation Intercalation\endlink).
 */
/*
 * Copyright (2005) Sandia Corporation. Under the terms of
 * Contract DE-AC04-94AL85000 with Sandia Corporation, the
 * U.S. Government retains certain rights in this software.
 */

//#include "Intercalation.h"
#include "cantera/thermo/Intercalation.h"
#include "cantera/thermo/PDSS.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/base/stringUtils.h"
#include "cantera/base/vec_functions.h"
#include "cantera/base/ctml.h"
#include "cantera/thermo/SpeciesThermoFactory.h"
#include "cantera/thermo/GeneralSpeciesThermo.h"

namespace Cantera
{
Intercalation::Intercalation(const Intercalation& right)
{
	*this = right;
}

Intercalation& Intercalation::operator=(const Intercalation& right)
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

ThermoPhase* Intercalation::duplMyselfAsThermoPhase() const
{
	return new Intercalation(*this);
}

int Intercalation::eosType() const
{
	return cIncompressible;
}

doublereal Intercalation::enthalpy_mole() const
{
	doublereal p0 = m_spthermo->refPressure();
	return GasConstant * temperature() *
			mean_X(enthalpy_RT()) + (pressure() - p0)/molarDensity();
}

doublereal Intercalation::entropy_mole() const
{
	return GasConstant * (mean_X(entropy_R()) - sum_xlogx());
}

doublereal Intercalation::cp_mole() const
{
	return GasConstant * mean_X(cp_R());
}

doublereal Intercalation::cv_mole() const
{
	return cp_mole();
}

doublereal Intercalation::pressure() const
{
	return m_press;
}

void Intercalation::setPressure(doublereal p)
{
	m_press = p;
}

void Intercalation::getActivityConcentrations(doublereal* c) const
{
	getConcentrations(c);
}

void Intercalation::getActivityCoefficients(doublereal* ac) const
{
	for (size_t k = 0; k < m_kk; k++) {
		ac[k] = 1.0;
	}
}
void Intercalation::getPartialMolarEnthalpies(doublereal* result) const {
	doublereal vdp = (pressure() - m_spthermo->refPressure()) / molarDensity();
	doublereal rt = temperature() * GasConstant;
	const vector_fp& h_RT = enthalpy_RT();
	for (size_t k = 0; k < m_kk; k++) {
		result[k] = rt*h_RT[k] + vdp;
	}
}

void Intercalation::getPartialMolarEntropies(doublereal* result) const {
	doublereal xx;
	doublereal r = GasConstant;
	const vector_fp& s_R = entropy_R();
	for (size_t k = 0; k < m_kk; k++) {
		xx = moleFraction(k);
		if(xx < SmallNumber) xx = SmallNumber;
		result[k] = r*(s_R[k] - log(xx));
	}
}
doublereal Intercalation::standardConcentration(size_t k) const
{
	return molarDensity();
}

void Intercalation::getChemPotentials(doublereal* mu) const
{
	doublereal vdp = (pressure() - m_spthermo->refPressure())/
			molarDensity();
	doublereal rt = temperature() * GasConstant;
	const vector_fp& g_RT = gibbs_RT();
	for (size_t k = 0; k < m_kk; k++) {
		double xx = std::max(SmallNumber, moleFraction(k));
		mu[k] = rt*(g_RT[k] + log(xx)) + vdp;
	}
}


void Intercalation::getStandardChemPotentials(doublereal* mu0) const
{
	getPureGibbs(mu0);
}

void Intercalation::initThermo()
{
	ThermoPhase::initThermo();
	m_h0_RT.resize(m_kk);
	m_g0_RT.resize(m_kk);
	m_cp0_R.resize(m_kk);
	m_s0_R.resize(m_kk);
	m_pp.resize(m_kk);
}

void Intercalation::setMoleFractions(const doublereal* const x)
{
	Phase::setMoleFractions(x);
	_updateThermo();
}

void Intercalation::setMoleFractions_NoNorm(const doublereal* const x)
{
	Phase::setMoleFractions_NoNorm(x);
	_updateThermo();
}

void Intercalation::setMassFractions(const doublereal* const y)
{
	Phase::setMassFractions(y);
	_updateThermo();
}

void Intercalation::setMassFractions_NoNorm(const doublereal* const y)
{
	Phase::setMassFractions_NoNorm(y);
	_updateThermo();
}

void Intercalation::setConcentrations(const doublereal* const conc)
{
	Phase::setConcentrations(conc);
	_updateThermo();
}

void Intercalation::setToEquilState(const doublereal* lambda_RT)
{
	throw CanteraError("setToEquilState","not yet impl.");
}

void Intercalation::_updateThermo() const
{
	doublereal tnow = temperature();
	doublereal xnow = moleFraction(m_kk_int);
	doublereal c[4];
	c[0]=tnow;
	c[1]=interp_h(xnow) * 1e3; // 1e3 for conversion J/mol -> J/kmol
	c[2]=interp_s(xnow) * 1e3 + GasConstant*std::log(xnow/(1.0-xnow)); // 1e3 for conversion J/K/mol -> J/K/kmol
	c[3]=0.0;

	if (m_tlast != tnow || m_xlast != xnow) {
		GeneralSpeciesThermo* genSpthermo = dynamic_cast<GeneralSpeciesThermo*>(m_spthermo);
		if (!genSpthermo) {
			throw CanteraError("Intercalation::_updateThermo",
					"unable to modify intercalation thermo");
		} else {
			genSpthermo->modifyParams(m_kk_int, &c[0]);
			m_spthermo->update(tnow, &m_cp0_R[0], &m_h0_RT[0],
					&m_s0_R[0]);
			for (size_t k = 0; k < m_kk; k++) {
				m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
			}
			m_xlast = xnow;
		}
		m_spthermo->update(tnow, &m_cp0_R[0], &m_h0_RT[0],
				&m_s0_R[0]);
		for (size_t k = 0; k < m_kk; k++) {
			m_g0_RT[k] = m_h0_RT[k] - m_s0_R[k];
		}
		m_tlast = tnow;
	}
}

void Intercalation::initThermoXML(XML_Node& phaseNode, const std::string& id_)
{
	Intercalation::initThermo();
	size_t offset;
	doublereal a=0,b=0,c=0;
	std::ifstream t;
	std::string line;

	if (phaseNode.hasChild("thermo")) {
		XML_Node& thermoNode = phaseNode.child("thermo");
		std::string model = thermoNode["model"];
		if (model == "Intercalation") {
			if (thermoNode.hasChild("data")) {
				XML_Node& scNode = thermoNode.child("data");
				std::vector<std::string> nameFile;
				getStringArray(scNode, nameFile);
				if (nameFile.size() != 1) {
					throw CanteraError("Intercalation::initThermoXML",
							"badly formed data XML node");
				}
				m_dataFile = nameFile[0];
			} else {
				throw CanteraError("initThermoXML",
						"Unspecified data file");
			}
			if (thermoNode.hasChild("intercalation_species")) {
				std::string intercalation_species_name = thermoNode.child("intercalation_species").value();
				m_kk_int = static_cast<size_t>(speciesIndex(intercalation_species_name));
				if (m_kk_int < 0) {
					throw CanteraError("Intercalation::initThermoXML",
							"Species " + intercalation_species_name + " not found.");
				}
			} else {
				throw CanteraError("Intercalation::initThermoXML",
						"Unspecified intercalation species");
			}
		} else {
			throw CanteraError("Intercalation::initThermoXML",
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
			std::sscanf(line.c_str(), "%lf %lf %lf", &a, &b, &c);
			molefrac_h.push_back(std::make_pair(a,b));
			molefrac_s.push_back(std::make_pair(a,c));
		}
	}
	t.close();
	ThermoPhase::initThermoXML(phaseNode, id_);
}

doublereal Intercalation::interp_h(doublereal x) const
{
	// Assumes that "table" is sorted by .first
	// Check if x is out of bound
	if (x > molefrac_h.back().first) return molefrac_h.back().second;
	if (x < molefrac_h[0].first) return molefrac_h[0].second;
	std::vector<std::pair<double, double> >::const_iterator it, it2;
	// INFINITY is defined in math.h in the glibc implementation
	it = lower_bound(molefrac_h.begin(), molefrac_h.end(), std::make_pair(x, molefrac_h[0].second));
	// Corner case
	if (it == molefrac_h.begin()) return it->second;
	it2 = it;
	--it2;
	return it2->second + (it->second - it2->second)*(x - it2->first)/(it->first - it2->first);
}

doublereal Intercalation::interp_s(doublereal x) const
{
	// Assumes that "table" is sorted by .first
	// Check if x is out of bound
	if (x > molefrac_s.back().first) return molefrac_s.back().second;
	if (x < molefrac_s[0].first) return molefrac_s[0].second;
	std::vector<std::pair<double, double> >::const_iterator it, it2;
	// INFINITY is defined in math.h in the glibc implementation
	it = lower_bound(molefrac_s.begin(), molefrac_s.end(), std::make_pair(x, molefrac_s[0].second));
	// Corner case
	if (it == molefrac_s.begin()) return it->second;
	it2 = it;
	--it2;
	return it2->second + (it->second - it2->second)*(x - it2->first)/(it->first - it2->first);
}

void Intercalation::setParametersFromXML(const XML_Node& eosdata)
{
	eosdata._require("model","Intercalation");
	doublereal rho = getFloat(eosdata, "density", "toSI");
	setDensity(rho);
}

}
