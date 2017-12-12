#include "gtest/gtest.h"
#include "cantera/thermo/IdealSolidSolnPhaseTabulatedThermo.h"
#include "cantera/thermo/ThermoFactory.h"

namespace Cantera
{

class IdealSolidSolnPhaseTabulatedThermo_Test : public testing::Test
{
public:
	IdealSolidSolnPhaseTabulatedThermo_Test(){
		test_phase.reset(newPhase("../data/IdealSolidSolnPhaseTabulatedThermo.cti"));
	}

	//vary the composition of a co2-h2 mixture:
    void set_defect_X(const double x) {
        vector_fp moleFracs(2);
        moleFracs[0] = x;
        moleFracs[1] = 1-x;
        test_phase->setMoleFractions(&moleFracs[0]);
    }
	
	std::unique_ptr<ThermoPhase> test_phase;
};

TEST_F(IdealSolidSolnPhaseTabulatedThermo_Test,construct_from_cti)
{
	IdealSolidSolnPhaseTabulatedThermo* IdealSolidSolnPhaseTabulatedThermo_phase = dynamic_cast<IdealSolidSolnPhaseTabulatedThermo*>(test_phase.get());
	EXPECT_TRUE(IdealSolidSolnPhaseTabulatedThermo_phase != NULL);
}

TEST_F(IdealSolidSolnPhaseTabulatedThermo_Test,interp_h)
{

	test_phase->setState_TP(298.15, 101325.);
	const double expected_result[9] = {
       -1915683.750000,
       -3140466.547093,
	   -4499406.744835,
	   -5758701.080648,
	   -6622802.262523,
	   -6028010.441154,
	   -6461378.467967,
	   -7198005.010302,
	   -7728562.111677
    };

    double xmin = 0.10;
    double xmax = 0.75;
    int numSteps= 9;
    double dx = (xmax-xmin)/(numSteps-1);
    for(int i=0; i < 9; ++i)
    {
        set_defect_X(xmin + i*dx);
        EXPECT_NEAR(expected_result[i], test_phase->enthalpy_mole(), 1.e-6);
    }

}

TEST_F(IdealSolidSolnPhaseTabulatedThermo_Test,interp_s)
{

	test_phase->setState_TP(298.15, 101325.);
	const double expected_result[9] = {
       	 965.047889,
	      30.260907,
	   -1835.910718,
	   -2110.486834,
	   -1871.952503,
	     745.480189,
	    1731.647583,
	    2773.311421,
	    4291.798734
    };

    double xmin = 0.10;
    double xmax = 0.75;
    int numSteps= 9;
    double dx = (xmax-xmin)/(numSteps-1);
    for(int i=0; i < 9; ++i)
    {
        set_defect_X(xmin + i*dx);
        EXPECT_NEAR(expected_result[i], test_phase->entropy_mole(), 1.e-6);
    }

}


TEST_F(IdealSolidSolnPhaseTabulatedThermo_Test,chem_potentials)
{
    test_phase->setState_TP(298.15,101325.);

    const double expected_result[9] = {
       -19683470.209968,
       -15137144.061807,
       -12934679.921407,
       -12928648.315684,
       -12413845.159368,
       -10639955.352516,
       -10335598.658919,
       -10643123.893675,
       -10865361.226623
    };

    double xmin = 0.10;
    double xmax = 0.75;
    int numSteps= 9;
    double dx = (xmax-xmin)/(numSteps-1);
    vector_fp chemPotentials(2);
    for(int i=0; i < 9; ++i)
    {
        set_defect_X(xmin + i*dx);
        test_phase->getChemPotentials(&chemPotentials[0]);
        EXPECT_NEAR(expected_result[i], chemPotentials[0], 1.e-6);
    }

}


TEST_F(IdealSolidSolnPhaseTabulatedThermo_Test,mole_fractions)
{
    test_phase->setState_TP(298.15,101325.);

    double xmin = 0.10;
    double xmax = 0.75;
    int numSteps= 9;
    double dx = (xmax-xmin)/(numSteps-1);
    vector_fp molefracs(2);
    for(int i=0; i < 9; ++i)
    {
        set_defect_X(xmin + i*dx);
        test_phase->getMoleFractions(&molefracs[0]);
        EXPECT_NEAR(xmin + i*dx, molefracs[0], 1.e-6);
    }

}

TEST_F(IdealSolidSolnPhaseTabulatedThermo_Test,partialMolarEntropies)
{
    test_phase->setState_TP(298.15,101325.);

    const double expected_result[9] = {
         1766.334764,
        -7343.846796,
       -14106.715418,
       -12825.546357,
       -10629.612570,
        -4250.296045,
        -2222.010231,
         -403.314051,
         1880.301004
    };

    double xmin = 0.10;
    double xmax = 0.75;
    int numSteps= 9;
    double dx = (xmax-xmin)/(numSteps-1);
    vector_fp partialMolarEntropies(2);
    for(int i=0; i < 9; ++i)
    {
        set_defect_X(xmin + i*dx);
        test_phase->getPartialMolarEntropies(&partialMolarEntropies[0]);
        EXPECT_NEAR(expected_result[i], partialMolarEntropies[0], 1.e-6);
    }

}


}