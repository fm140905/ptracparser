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
    MCNPptrac = new MCNPPTRACBinary("/media/ming/DATA/projects/TRISO/test3/ptracparser/data/slabbinp");
  }

  void TearDown()
  {
    delete MCNPptrac;
  }
};

TEST_F(MCNPtestPtracBinary, ReadFirstData)
{
  MCNPptrac->readNextPtracData(1000);
  auto const &record = MCNPptrac->getPTRACRecord();
  ASSERT_EQ(record.size(), 11);
  // ASSERT_EQ(record.pointID, 1);
  // ASSERT_EQ(record.eventID, 1000);
  // ASSERT_EQ(record.cellID, 3001);
  // ASSERT_DOUBLE_EQ(record.point[0], 12.02427913688436);
  // ASSERT_DOUBLE_EQ(record.point[1], -72.881828858643473);
  // ASSERT_DOUBLE_EQ(record.point[2], 1.0882843926953107);
}

// TEST_F(MCNPtestPtracBinary, ReadAll)
// {
//   int ii = 1;
//   while (ii <= 1000) {
//     MCNPptrac->readNextPtracData(2000);
//     ii++;
//   }

//   auto const &record = MCNPptrac->getPTRACRecord();
//   ASSERT_EQ(record.pointID, 1000);
//   ASSERT_EQ(record.eventID, 1000);

//   ASSERT_EQ(record.cellID, 1001);
//   ASSERT_EQ(record.materialID, 3);

//   ASSERT_DOUBLE_EQ(record.point[0], -2.2580291076456316);
//   ASSERT_DOUBLE_EQ(record.point[1], 18.880472384323852);
//   ASSERT_DOUBLE_EQ(record.point[2], -1.2855986393557863);

//   ASSERT_EQ(MCNPptrac->getNbPointsRead(), 1000);
// }


