#include "crpropa/PhotonBackground.h"
#include "crpropa/Candidate.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"

#include <stdlib.h>

#include "gtest/gtest.h"

#include <fstream>

namespace crpropa {

TEST(PhotonField, initNativeFields) {
	// test if native field files can be found and initialized
	std::string nativeFieldNames[10] = {"CMB",
										"IRB_Kneiske04",
										"IRB_Stecker05",
										"IRB_Finke10",
										"IRB_Dominguez11",
										"IRB_Franceschini08",
										"IRB_Gilmore12",
										"IRB_Stecker16_lower",
										"IRB_Stecker16_upper",
										"URB_Protheroe96"};
	for (int i = 0; i < 10; ++i) {
		std::string fieldLoc = "Scaling/" + std::string(nativeFieldNames[i]) + ".txt";
		CustomPhotonField cpf = CustomPhotonField(getDataPath(fieldLoc));
	}
}

TEST(PhotonField, nativeFieldFilesFormat) {
	// test if native field files have the correct format
	std::string nativeFieldNames[10] = {"CMB",
										"IRB_Kneiske04",
										"IRB_Stecker05",
										"IRB_Finke10",
										"IRB_Dominguez11",
										"IRB_Franceschini08",
										"IRB_Gilmore12",
										"IRB_Stecker16_lower",
										"IRB_Stecker16_upper",
										"URB_Protheroe96"};
	for (int i = 0; i < 10; ++i) {
		std::string fieldLoc = "Scaling/" + std::string(nativeFieldNames[i]) + ".txt";
		CustomPhotonField cpf = CustomPhotonField(getDataPath(fieldLoc));
		EXPECT_GT(cpf.photonEnergy.size(), 0);
		EXPECT_GT(cpf.photonRedshift.size(), 0);
		EXPECT_GT(cpf.photonDensity.size(), 0);
		EXPECT_EQ(cpf.photonEnergy.size() * cpf.photonRedshift.size(), cpf.photonDensity.size());
	}
}


TEST(PhotonField, customFieldFilesFormat) {
	// test if custom field files can be initialized correctly (if any)
	// assume that field slot is unused if initialization fails
	std::string customFieldNames[8] = {"PF1", "PF2", "PF3", "PF4",
									   "PF5", "PF6", "PF7", "PF8"};
	for (int i = 0; i < 8; ++i) {
		std::string fieldLoc = "Scaling/" + std::string(customFieldNames[i]) + ".txt";
		try
		{	
			CustomPhotonField cpf = CustomPhotonField(getDataPath(fieldLoc));
			EXPECT_GT(cpf.photonEnergy.size(), 0);
			EXPECT_GT(cpf.photonRedshift.size(), 0);
			EXPECT_GT(cpf.photonDensity.size(), 0);
			EXPECT_EQ(cpf.photonEnergy.size() * cpf.photonRedshift.size(), cpf.photonDensity.size());
		}
		catch (std::runtime_error &e)
		{
			continue;  // assume that field slot is unused if initialization fails
		}
	}		
}

TEST(PhotonField, customFieldFilesEntries) {
	// test if custom field files comply with the required format
	// assume that field slot is unused if initialization fails
	std::string customFieldNames[8] = {"PF1", "PF2", "PF3", "PF4",
									   "PF5", "PF6", "PF7", "PF8"};
	for (int i = 0; i < 8; ++i) {
		std::string fieldLoc = "Scaling/" + std::string(customFieldNames[i]) + ".txt";
		try
		{	
			CustomPhotonField cpf = CustomPhotonField(getDataPath(fieldLoc));
			// test if photon energies are in ascending order
			for (int i = 1; i < cpf.photonEnergy.size(); ++i) {
				EXPECT_LT(cpf.photonEnergy[i-1], cpf.photonEnergy[i]);
			}
			// test if photon redshifts are in ascending order
			for (int i = 1; i < cpf.photonRedshift.size(); ++i) {
				EXPECT_LT(cpf.photonRedshift[i-1], cpf.photonRedshift[i]);
			}
			// test if photon densities are >= 0
			for (int i = 0; i < cpf.photonDensity.size(); ++i) {
				EXPECT_GE(cpf.photonDensity[i], 0.);
			}
		}
		catch (std::runtime_error &e)
		{
			continue;  // assume that field slot is unused if initialization fails
		}
	}		
}

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}

} // namespace crpropa
