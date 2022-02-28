#include "pulse.hh"
#include "gtest/gtest.h"
#include <fstream>

class PulseTest : public ::testing::Test
{
public:
    MCNPPTRAC *MCNPptrac;
    void SetUp()
    {
        MCNPptrac = new MCNPPTRACBinary("../../tests/data/MCNP6.2/ptrac");
    }

    void TearDown()
    {
        delete MCNPptrac;
    }
};

TEST_F(PulseTest, FirstPulse)
{
    MCNPptrac->readNextNPS(1000);
    const NPSHistory& record = MCNPptrac->getNPSHistory();

    // first particle
    ParticleHistory parHist = record[0];
    // first pulse
    Pulse firstPulse = Pulse(parHist);
    Pulse truth;
    truth.startPos = std::vector<double>{9.76040, -0.81083, 9.02960};
    truth.endPos = std::vector<double>{ 9.63740, -0.75003, 8.76890};
    truth.energy = 0.31638;
    truth.time = 7.844400E+07;
    // truth.latticeIndex = std::make_pair(21, 33);
    truth.nps = 1;

    for (int i = 0; i < 3; i++)
    {
        EXPECT_NEAR(firstPulse.startPos[i], truth.startPos[i], 0.001);
    }
    for (int i = 0; i < 3; i++)
    {
        EXPECT_NEAR(firstPulse.endPos[i], truth.endPos[i], 0.001);
    }
    EXPECT_NEAR(firstPulse.energy, truth.energy, 0.001);
    EXPECT_NEAR(firstPulse.time, truth.time, truth.time * 0.001);
    // EXPECT_EQ(firstPulse.latticeIndex.first, truth.latticeIndex.first);
    // EXPECT_EQ(firstPulse.latticeIndex.second, truth.latticeIndex.second);
}

TEST_F(PulseTest, MidPulse)
{
    const int pulseNum = 100;
    int pulseIdx = 0;
    Pulse midPulse;
    while (MCNPptrac->readNextNPS(1000))
    {
        const NPSHistory record = MCNPptrac->getNPSHistory();
        for (auto iter = record.begin(); iter != record.end(); iter++)
        {
            midPulse = Pulse(*iter);
            if (midPulse.energy > 0)
            {
                // vaild pulse
                pulseIdx++;
                if (pulseIdx >= pulseNum)
                {
                    break;
                }
            }
        }
        if (pulseIdx >= pulseNum)
        {
            break;
        }
    }
    Pulse truth;
    truth.startPos = std::vector<double>{9.35390, -7.20160, 7.26640};
    truth.endPos = std::vector<double>{9.47970, -6.99500, 7.49860};
    truth.energy = 0.39276;
    truth.time = 3.612400E+07;
    // truth.latticeIndex = std::make_pair(13, 37);
    truth.nps = 105;

    for (int i = 0; i < 3; i++)
    {
        EXPECT_NEAR(midPulse.startPos[i], truth.startPos[i], 0.001);
    }
    for (int i = 0; i < 3; i++)
    {
        EXPECT_NEAR(midPulse.endPos[i], truth.endPos[i], 0.001);
    }
    EXPECT_NEAR(midPulse.energy, truth.energy, 0.001);
    EXPECT_NEAR(midPulse.time, truth.time, truth.time * 0.001);
    // EXPECT_EQ(midPulse.latticeIndex.first, truth.latticeIndex.first);
    // EXPECT_EQ(midPulse.latticeIndex.second, truth.latticeIndex.second);
}


TEST_F(PulseTest, LastPulse)
{
    const int pulseNum = 1e6;
    int pulseIdx = 0;
    Pulse lastPulse;
    while (MCNPptrac->readNextNPS(1e6))
    {
        const NPSHistory record = MCNPptrac->getNPSHistory();
        for (auto iter = record.begin(); iter != record.end(); iter++)
        {
            lastPulse = Pulse(*iter);
            if (lastPulse.energy > 0)
            {
                // vaild pulse
                pulseIdx++;
                if (pulseIdx >= pulseNum)
                {
                    break;
                }
            }
        }
        if (pulseIdx >= pulseNum)
        {
            break;
        }
    }
    Pulse truth;
    truth.startPos = std::vector<double>{-9.93780, 4.93200, 2.04840};
    truth.endPos = std::vector<double>{-10.21900, 4.76060, 2.08380};
    truth.energy = 0.50319;
    truth.time = 7.608400E+06;
    // truth.latticeIndex = std::make_pair(28, 8);
    truth.nps = 45609;

    for (int i = 0; i < 3; i++)
    {
        EXPECT_NEAR(lastPulse.startPos[i], truth.startPos[i], 0.001);
    }
    for (int i = 0; i < 3; i++)
    {
        EXPECT_NEAR(lastPulse.endPos[i], truth.endPos[i], 0.001);
    }
    EXPECT_NEAR(lastPulse.energy, truth.energy, 0.001);
    EXPECT_NEAR(lastPulse.time, truth.time, truth.time * 0.001);
    // EXPECT_EQ(lastPulse.latticeIndex.first, truth.latticeIndex.first);
    // EXPECT_EQ(lastPulse.latticeIndex.second, truth.latticeIndex.second);
}