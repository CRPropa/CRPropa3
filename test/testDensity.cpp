#include "crpropa/massDistribution/Massdistribution.h"
#include "crpropa/massDistribution/Cordes.h"
#include "crpropa/massDistribution/Ferriere.h"
#include "crpropa/massDistribution/Nakanishi.h"
#include "crpropa/massDistribution/ConstantDensity.h"
#include "crpropa/Units.h"

#include "gtest/gtest.h"
#include <stdexcept>
#include <cmath>

namespace crpropa {

TEST(testConstantDensity, SimpleTest) {
	//test ConstantDensity in all types and in total density (output)
	ConstantDensity n(2/ccm,3/ccm, 2/ccm);
	Vector3d p(1*pc,2*pc,1*kpc); 	// random position for testing density
	EXPECT_DOUBLE_EQ(n.getHIDensity(p), 2e6);	// density output in m^-3
	EXPECT_DOUBLE_EQ(n.getHIIDensity(p), 3e6);
	EXPECT_DOUBLE_EQ( n.getH2Density(p), 2e6);
	EXPECT_DOUBLE_EQ(n.getDensity(p), 7e6);		// total density 2+3+2 = 7 (/ccm)
	EXPECT_DOUBLE_EQ(n.getNucleonDensity(p),9e6);	// nucleon density 2+3+2*2 = 9 (/ccm) factor 2 for molecular hydrogen
	
	//test set/get function for used type
	bool useHI = n.getIsForHI();
	bool useHII= n.getIsForHII();
	bool useH2 = n.getIsForH2();
	//check if all types are activited
	EXPECT_TRUE(useHI); 
	EXPECT_TRUE(useHII);
	EXPECT_TRUE(useH2);
	
	//set density number to 500
	n.setHI(500.);
	n.setHII(500.);
	n.setH2(500.);
		
	//check if output is changed to 500 in types and is 0 (all deactivated) for Hges
	EXPECT_DOUBLE_EQ(n.getHIDensity(p), 500.);	
	EXPECT_DOUBLE_EQ(n.getHIIDensity(p), 500.);
	EXPECT_DOUBLE_EQ(n.getH2Density(p), 500.);
	
	//deactivate all 
	n.setHI(false);
	n.setHII(false);
	n.setH2(false);
	
	//check if all isfor are set to false
	useHI = n.getIsForHI();
	useHII= n.getIsForHII();
	useH2 = n.getIsForH2();
	EXPECT_FALSE(useHI);
	EXPECT_FALSE(useHII);
	EXPECT_FALSE(useH2);
	
	//get<type>Density is independent of type activation, getDensity is not independent 
	//check if getDensity returns 0. (should give a error message in log file)
	EXPECT_DOUBLE_EQ(n.getDensity(p), 0);
	EXPECT_DOUBLE_EQ(n.getNucleonDensity(p),0);	
}

TEST(testCustomDensity, SimpleTest) {
	
	CustomDensity MD;
	
	Vector3d p(50*pc, 20*pc, -100*pc);	//random position for testing density 
	
	//try to get density without load any option in. Should give Warning in log-File and return density of 0.
	EXPECT_DOUBLE_EQ(MD.getDensity(p),0);	
	
	MD.add(new ConstantDensity(1,1,1));	//all types loaded
	MD.add(new ConstantDensity(0,2,0));	//overwrite HII type 
	MD.add(new ConstantDensity(0,0,3));	//overwrite H2 type
	
	//check get density output
	EXPECT_DOUBLE_EQ(MD.getHIDensity(p),1);
	EXPECT_DOUBLE_EQ(MD.getHIIDensity(p),2);
	EXPECT_DOUBLE_EQ(MD.getH2Density(p),3);
	EXPECT_DOUBLE_EQ(MD.getDensity(p),6);	//total density 2+2+3
	EXPECT_DOUBLE_EQ(MD.getNucleonDensity(p),9);	// nucleon density = 2 + 2 + 2*3  factor 2 for molecular hydrogen
	
	//check get function for type 
	EXPECT_TRUE(MD.getIsForHI());
	EXPECT_TRUE(MD.getIsForHII());
	EXPECT_TRUE(MD.getIsForH2());
	
	//check type deactivate funktion
	MD.setIsForHI(false);
	EXPECT_FALSE(MD.getIsForHI());
	EXPECT_TRUE(MD.getIsForHII());
	EXPECT_TRUE(MD.getIsForH2());
	
	MD.setIsForHII(false);
	EXPECT_FALSE(MD.getIsForHI());
	EXPECT_FALSE(MD.getIsForHII());
	EXPECT_TRUE(MD.getIsForH2());
	
	MD.setIsForH2(false);
	EXPECT_FALSE(MD.getIsForHI());
	EXPECT_FALSE(MD.getIsForHII());
	EXPECT_FALSE(MD.getIsForH2());
	
	//check density output if all types are deactivated
	//should give error message in log-file and return density of 0.
	EXPECT_DOUBLE_EQ(MD.getDensity(p),0.);
	EXPECT_DOUBLE_EQ(MD.getNucleonDensity(p),0.);

} 

TEST(testDensityList, SimpleTest) {

	DensityList MS;
	MS.addDensity(new ConstantDensity(1,1,2));	//sum 4
	MS.addDensity(new ConstantDensity(2,3,1));	//sum 6
	
	Vector3d p(50*pc,10*pc,-30*pc);	//random position for testing density
	EXPECT_DOUBLE_EQ(MS.getHIDensity(p),3);
	EXPECT_DOUBLE_EQ(MS.getHIIDensity(p),4);
	EXPECT_DOUBLE_EQ(MS.getH2Density(p),3);
	EXPECT_DOUBLE_EQ(MS.getDensity(p),10);	//sum of sums
	EXPECT_DOUBLE_EQ(MS.getNucleonDensity(p),13); 	// 3+4+2*3 factor 2 for molecular hydrogen
	
}

TEST(testCordes, checkValueAtCertainPoints) {

	Cordes n;
	
	//check type Information
	EXPECT_FALSE(n.getIsForHI());
	EXPECT_TRUE(n.getIsForHII());
	EXPECT_FALSE(n.getIsForH2());

	Vector3d p(3.1*kpc,2.9*kpc,-30*pc);	//position for testing density	
	
	EXPECT_NEAR(n.getHIIDensity(p), 184500.,1);	// output in m^-3 ; uncertainty of 1e-6 cm^-3 
	EXPECT_NEAR(n.getDensity(p), 184500.,1);	
	EXPECT_NEAR(n.getNucleonDensity(p), 184500,1);	// only HII component -> no differenz between density and nucleon density
	p.z=30*pc;			// invariant density for +/- z
	EXPECT_NEAR(n.getDensity(p),184500,1); 
	
}

TEST(testNakanishi, checkValueAtCertainPoints) {

	Nakanishi n;
	
	//check type Information
	EXPECT_TRUE(n.getIsForHI());
	EXPECT_FALSE(n.getIsForHII());
	EXPECT_TRUE(n.getIsForH2());
	
	//first position for testing density
	Vector3d p(4*kpc,-2.5*kpc,-0.85*kpc);	
	
	//testing HI component
	EXPECT_NEAR(n.getHIPlanedensity(p),162597,1);	//uncertaincy of 1e-6 cm^-3 
	EXPECT_NEAR(n.getHIScaleheight(p),0.3109*kpc,0.1*pc); 
	EXPECT_NEAR(n.getHIDensity(p),914,1); 	//uncertainc 1e-6 cm^-3 
	
	//testing H2 compontent
	EXPECT_NEAR(n.getH2Planedensity(p),741999,1); //uncertaincy of 1e-6 cm^-3
	EXPECT_NEAR(n.getH2Scaleheight(p),88.2*pc,0.1*pc);
	EXPECT_NEAR(n.getH2Density(p),0,1);
	
	//testing total Density
	EXPECT_NEAR(n.getDensity(p),914,2); //double uncertaincy for both type á 1cm^-3
	EXPECT_NEAR(n.getNucleonDensity(p),914,2);	// 914 + 0*2	factor 2 for molecular hydrogen	
	
	//second position for testing density
	p = Vector3d(50*pc,100*pc,10*pc);	
	
	//testing HI component
	EXPECT_NEAR(n.getHIPlanedensity(p),543249,1);
	EXPECT_NEAR(n.getHIScaleheight(p),125.6*pc,0.1*pc);
	EXPECT_NEAR(n.getHIDensity(p),540867,1);
	
	//testing H2 component
	EXPECT_NEAR(n.getH2Planedensity(p),10556748,1);
	EXPECT_NEAR(n.getH2Scaleheight(p),57.2*pc,0.1*pc);
	EXPECT_NEAR(n.getH2Density(p),10335137,1);
	
	//testing total Density
	EXPECT_NEAR(n.getDensity(p),10876004,2);
	EXPECT_NEAR(n.getNucleonDensity(p),21211141,2);	// factor 2 in molecular hydrogen
	
	//test set type function
	n.setIsForHI(false);
	EXPECT_FALSE(n.getIsForHI());
	EXPECT_TRUE(n.getIsForH2());
	
	n.setIsForH2(false);
	EXPECT_FALSE(n.getIsForHI());
	EXPECT_FALSE(n.getIsForH2());
	
	//check if density output is zero if all density-types are deaktivated (should give warning in log-file)
	EXPECT_DOUBLE_EQ(n.getDensity(p),0);	
	EXPECT_DOUBLE_EQ(n.getNucleonDensity(p),0);
	
}

TEST(testFerriere, checkValueAtCertainPoints) {
	
	Ferriere n;
	
	//check type information
	
	EXPECT_TRUE(n.getIsForHI());
	EXPECT_TRUE(n.getIsForHII());
	EXPECT_TRUE(n.getIsForH2());
	
	//testing density in inner Ring (R <= 3*kpc)
	Vector3d p(60*pc,-60*pc,-20*pc);	//testing position in region of CMZ
	
	//test CMZ Trafo
	Vector3d Trafo;
	Trafo = n.CMZTrafo(p);
	EXPECT_NEAR(Trafo.x,5.9767*pc,1e-4*pc);
	EXPECT_NEAR(Trafo.y,12.8171*pc,1e-4*pc);
	EXPECT_DOUBLE_EQ(Trafo.z,p.z);	//no transformation in z component
	
	//test DISK Trafo
	Trafo = n.DISKTrafo(p);
	EXPECT_NEAR(Trafo.x,11.0660*pc,1e-4*pc);
	EXPECT_NEAR(Trafo.y,82.5860*pc,1e-4*pc);
	EXPECT_NEAR(Trafo.z,-25.6338*pc,1e-4*pc);
	
	//testing density 
	EXPECT_NEAR(n.getHIDensity(p),6237723,1); 	//uncertaincy 1e-6 cm^-3
	EXPECT_NEAR(n.getH2Density(p),35484825,1); 
	EXPECT_NEAR(n.getHIIDensity(p),6243793,1);	
	EXPECT_NEAR(n.getDensity(p),47966341,1);
	EXPECT_NEAR(n.getNucleonDensity(p),83451166,2);		//factor 2 in molecular hydrogen; double uncertaincy
	
	Vector3d p2(-500*pc,-900*pc,35*pc);	//testing position in region of the DISK
	EXPECT_NEAR(n.getHIIDensity(p2),48190,1);
	EXPECT_NEAR(n.getHIDensity(p2),5,1);
	EXPECT_NEAR(n.getH2Density(p2),0,1);
	EXPECT_NEAR(n.getDensity(p2),48195,1);
	EXPECT_NEAR(n.getNucleonDensity(p2),48195,1);	// no H2 component -> no difference between density and nucleon-density
	
	//testing the outer region R>3kpc
	Vector3d p3(5*kpc,4*kpc,-29*pc);	//testing position with 3kpc < R < R_sun
	EXPECT_NEAR(n.getHIDensity(p3),540607,1);
	EXPECT_NEAR(n.getHIIDensity(p3),66495 ,1);
	EXPECT_NEAR(n.getH2Density(p3),2492685,1);
	EXPECT_NEAR(n.getDensity(p3),3099787,1);
	EXPECT_NEAR(n.getNucleonDensity(p3), 5592472,1);	
	
	Vector3d p4(10*kpc,2*kpc,50*pc);	//testing position with R > R_sun
	EXPECT_NEAR(n.getHIDensity(p4),431294,1);
	EXPECT_NEAR(n.getHIIDensity(p4),22109,1);
	EXPECT_NEAR(n.getH2Density(p4),54099,1);
	EXPECT_NEAR(n.getDensity(p4),507502,1);
	EXPECT_NEAR(n.getNucleonDensity(p4),561601,1);
	
	//test get/set type function
	n.setIsForHI(false);
	EXPECT_FALSE(n.getIsForHI());
	EXPECT_TRUE(n.getIsForHII());
	EXPECT_TRUE(n.getIsForH2());
	
	n.setIsForHII(false);
	EXPECT_FALSE(n.getIsForHI());
	EXPECT_FALSE(n.getIsForHII());
	EXPECT_TRUE(n.getIsForH2());
	
	n.setIsForH2(false);
	EXPECT_FALSE(n.getIsForHI());
	EXPECT_FALSE(n.getIsForHII());
	EXPECT_FALSE(n.getIsForH2());
	
	//check if density is set to zero if all types are deactivated (should give warning in log-file)
	EXPECT_DOUBLE_EQ(n.getDensity(p),0);
	EXPECT_DOUBLE_EQ(n.getNucleonDensity(p),0);
}
} //namespace crpropa
