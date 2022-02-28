/**
* @file parser_test.cc
*
*
* @brief Test of binary PTRAC file parser
*
* @author jofausti, Ming Fang
* @version 1.1
*/
#include "parser.hh"
#include "gtest/gtest.h"
#include <fstream>

using namespace std;

class MCNPtestPtracBinary : public ::testing::Test
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

TEST_F(MCNPtestPtracBinary, ReadFirstData)
{
  MCNPptrac->readNextNPS(1000);
  const NPSHistory& record = MCNPptrac->getNPSHistory();
  EXPECT_EQ(record.size(), 4);

  // first particle
  ParticleHistory parHist = record[0];
  EXPECT_EQ(parHist.size(), 3);
  // first event
  const auto beg = parHist.begin();
  const Event firstEvent = Event({1, 2030, 602, 
                              {0.97606E+01, -0.81088E+00, 0.90298E+01}, 
                              0.14663E+01, 0.10000E+01, 0.78444E+08});
  EXPECT_EQ(beg->nps, firstEvent.nps);
  EXPECT_EQ(beg->eventID, firstEvent.eventID);
  EXPECT_EQ(beg->cellID, firstEvent.cellID);
  for (int i = 0; i < 3; i++)
  {
    EXPECT_NEAR(beg->pos[i], firstEvent.pos[i], 0.001);
  }
  EXPECT_NEAR(beg->energy, firstEvent.energy, 0.001);
  EXPECT_NEAR(beg->weight, firstEvent.weight, 0.001);
  EXPECT_NEAR(beg->time, firstEvent.time, 0.001 * beg->time);

  // last particle
  parHist = record[record.size() - 1];
  EXPECT_EQ(parHist.size(), 3);
  // last event
  const auto endd = parHist.rbegin();
  const Event lastEvent = Event({1, 5000, 603, 
                              {-0.68334E+01, 0.11587E+02, 0.55330E+01}, 
                              0.10000E-02, 0.10000E+01, 0.78441E+08});
  EXPECT_EQ(endd->nps, lastEvent.nps);
  EXPECT_EQ(endd->eventID, lastEvent.eventID);
  EXPECT_EQ(endd->cellID, lastEvent.cellID);
  for (int i = 0; i < 3; i++)
  {
    EXPECT_NEAR(endd->pos[i], lastEvent.pos[i], 0.001);
  }
  EXPECT_NEAR(endd->energy, lastEvent.energy, 0.001);
  EXPECT_NEAR(endd->weight, lastEvent.weight, 0.001);
  EXPECT_NEAR(endd->time, lastEvent.time, 0.001 * beg->time);

}

TEST_F(MCNPtestPtracBinary, ReadMiddle)
{
  int ii = 0;
  while (ii < 100) {
    MCNPptrac->readNextNPS(100);
    ii++;
  }

  // the 100-th nps
  const NPSHistory& record = MCNPptrac->getNPSHistory();
  EXPECT_EQ(record.size(), 3);
  
  // first particle
  ParticleHistory parHist = record[0];
  EXPECT_EQ(parHist.size(), 4);
  // first event
  const auto beg = parHist.begin();
  const Event firstEvent = Event({198, 2030, 602, 
                              {0.12243E+01, 0.11980E+02, -0.16554E+01}, 
                              0.14663E+01, 0.10000E+01, 0.26626E+08});
  EXPECT_EQ(beg->nps, firstEvent.nps);
  EXPECT_EQ(beg->eventID, firstEvent.eventID);
  EXPECT_EQ(beg->cellID, firstEvent.cellID);
  for (int i = 0; i < 3; i++)
  {
    EXPECT_NEAR(beg->pos[i], firstEvent.pos[i], 0.001);
  }
  EXPECT_NEAR(beg->energy, firstEvent.energy, 0.001);
  EXPECT_NEAR(beg->weight, firstEvent.weight, 0.001);
  EXPECT_NEAR(beg->time, firstEvent.time, 0.001 * beg->time);

  // last particle
  parHist = record[record.size() - 1];
  EXPECT_EQ(parHist.size(), 4);
  // last event
  const auto endd = parHist.rbegin();
  const Event lastEvent = Event({198, 5000, 602, 
                              {0.14093E+02, 0.37142E+01, 0.57761E+01}, 
                              0.10000E-02, 0.10000E+01, 0.26622E+08});
  EXPECT_EQ(endd->nps, lastEvent.nps);
  EXPECT_EQ(endd->eventID, lastEvent.eventID);
  EXPECT_EQ(endd->cellID, lastEvent.cellID);
  for (int i = 0; i < 3; i++)
  {
    EXPECT_NEAR(endd->pos[i], lastEvent.pos[i], 0.001);
  }
  EXPECT_NEAR(endd->energy, lastEvent.energy, 0.001);
  EXPECT_NEAR(endd->weight, lastEvent.weight, 0.001);
  EXPECT_NEAR(endd->time, lastEvent.time, 0.001 * beg->time);
}

TEST_F(MCNPtestPtracBinary, ReadLast)
{
  while (MCNPptrac->readNextNPS(1e6)) {

  }

  // the last nps
  const NPSHistory& record = MCNPptrac->getNPSHistory();
  EXPECT_EQ(record.size(), 3);
  
  // first particle
  ParticleHistory parHist = record[0];
  EXPECT_EQ(parHist.size(), 2);
  // first event
  const auto beg = parHist.begin();
  const Event firstEvent = Event({45609, 2030, 602, 
                              {-0.66029E+01, -0.86944E+01, 0.10149E+01}, 
                              0.14663E+01, 0.10000E+01, 0.76030E+07});
  EXPECT_EQ(beg->nps, firstEvent.nps);
  EXPECT_EQ(beg->eventID, firstEvent.eventID);
  EXPECT_EQ(beg->cellID, firstEvent.cellID);
  for (int i = 0; i < 3; i++)
  {
    EXPECT_NEAR(beg->pos[i], firstEvent.pos[i], 0.001);
  }
  EXPECT_NEAR(beg->energy, firstEvent.energy, 0.001);
  EXPECT_NEAR(beg->weight, firstEvent.weight, 0.001);
  EXPECT_NEAR(beg->time, firstEvent.time, 0.001 * beg->time);

  // last particle
  parHist = record[record.size() - 1];
  EXPECT_EQ(parHist.size(), 4);
  // last event
  const auto endd = parHist.rbegin();
  const Event lastEvent = Event({45609, 5000, 602, 
                              {-0.10219E+02, 0.47605E+01, 0.20838E+01}, 
                              0.10000E-02, 0.10000E+01, 0.76084E+07});
  EXPECT_EQ(endd->nps, lastEvent.nps);
  EXPECT_EQ(endd->eventID, lastEvent.eventID);
  EXPECT_EQ(endd->cellID, lastEvent.cellID);
  for (int i = 0; i < 3; i++)
  {
    EXPECT_NEAR(endd->pos[i], lastEvent.pos[i], 0.001);
  }
  EXPECT_NEAR(endd->energy, lastEvent.energy, 0.001);
  EXPECT_NEAR(endd->weight, lastEvent.weight, 0.001);
  EXPECT_NEAR(endd->time, lastEvent.time, 0.001 * beg->time);
}


