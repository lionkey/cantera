#include "gtest/gtest.h"
#include "cantera/thermo/ConstDensityTabulatedThermo.h"
#include "cantera/thermo/ThermoFactory.h"

namespace Cantera
{

class ConstDensityTabulatedThermo_Test : public testing::Test
{
public:
	ConstDensityTabulatedThermo_Test(){
		test_phase.reset(newPhase("../data/ConstDensityTabulatedThermo.cti"));
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

TEST_F(ConstDensityTabulatedThermo_Test,construct_from_cti)
{
	ConstDensityTabulatedThermo* constDensityTabulatedThermo_phase = dynamic_cast<ConstDensityTabulatedThermo*>(test_phase.get());
	EXPECT_TRUE(constDensityTabulatedThermo_phase != NULL);
}

TEST_F(ConstDensityTabulatedThermo_Test,interp_h)
{

	test_phase->setState_TP(298.15, 101325.);
	const double expected_result[9] = {
       	-1523717.129950,
       	-2332251.932434,
	    -3630175.923912,
	    -4707368.353702,
	    -5177756.746737,
	    -4674479.634792,
	    -5341855.069058,
	    -5661616.443350,
	    -4929484.137938
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

TEST_F(ConstDensityTabulatedThermo_Test,interp_s)
{

	test_phase->setState_TP(298.15, 101325.);
	const double expected_result[9] = {
       	1612.096614,
	    1065.699430,
	    -575.958400,
	   -1188.173091,
	    -873.828098,
	    2286.427038,
	    3388.709448,
	    4285.114152,
	    3749.463175
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

}