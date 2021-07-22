#include "crpropa/DiffusionTensor.h"

#include "gtest/gtest.h"

using namespace crpropa;

TEST(testQLTDiffusion, TestGetSetFunction){
    QLTDiffusion* Tens = new QLTDiffusion();

    // check default values 
    EXPECT_DOUBLE_EQ(Tens-> getEpsilon(), 0.1);
    EXPECT_DOUBLE_EQ(Tens-> getKappa0(), 6.1e24);
    EXPECT_DOUBLE_EQ(Tens-> getAlpha(), 1./3.);

    // set new values and check them
    Tens -> setEpsilon(0.3);
    EXPECT_DOUBLE_EQ(Tens-> getEpsilon(), 0.3);
    Tens -> setKappa0(1e5);
    EXPECT_DOUBLE_EQ(Tens-> getKappa0(), 1e5);
    Tens -> setAlpha(0.2);
    EXPECT_DOUBLE_EQ(Tens-> getAlpha(), 0.2);
}

TEST(testQLTDiffusion, TestForEnergyScaling){
    QLTDiffusion* Tens = new QLTDiffusion();

    int id = 1000010010; // protons

    // checking at norm Value 4 GeV
    double energy = 4*GeV;
    Candidate cand(id, energy);
    double norm = Tens-> getKappa0();
    EXPECT_DOUBLE_EQ(Tens-> getKappaParallel(&cand), norm);
    EXPECT_DOUBLE_EQ(Tens-> getKappaPerpendicular(&cand), 0.1*norm); // with epsilon = 0.1 as default
    EXPECT_DOUBLE_EQ(Tens-> getKappaPerpendicular2(&cand), 0.1*norm);

    // checking for another energy
    energy = 15*GeV;
    cand.current.setEnergy(energy);
    EXPECT_NEAR(Tens-> getKappaParallel(&cand), 9.47706e24, 1e19);
    EXPECT_NEAR(Tens-> getKappaPerpendicular(&cand), 9.47706e23, 1e18);
    EXPECT_NEAR(Tens-> getKappaPerpendicular2(&cand), 9.47706e23, 1e18);

    // checking energy far away from norming
    energy = 15*PeV;
    cand.current.setEnergy(energy);
    EXPECT_NEAR(Tens-> getKappaParallel(&cand), 9.47706e26, 1e21);
    EXPECT_NEAR(Tens-> getKappaPerpendicular(&cand),  9.47706e25, 1e20);
    EXPECT_NEAR(Tens-> getKappaPerpendicular2(&cand),  9.47706e25, 1e20);
}

