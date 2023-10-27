#pragma once

#pragma once


#include <array>

#include "st_tree/st_tree.h"


#define IF_REACTION 0
#define SF_REACTION 1
#define N2N_REACTION 2
#define B10_CAPTURE 3
#define H1_CAPTURE 4
#define CF252_CAPTURE 5
#define N3N_REACTION 6
#define N4N_REACTION 7

struct NeutronHistory {
  long nps;
  int creation_reaction; // 0 for IF, 1 for SF, 2 for (n,2n)
  // long creation_event_ID;
  std::array<double,3> creation_pos;
  double creation_energy;
  double creation_time;

  int destruction_reaction; // 0 for IF, 1 for SF, 2 for (n,2n), 
  // long destruction_event_ID;
  std::array<double,3> destruction_pos;
  double destruction_energy;
  double destruction_time;
  // 3 for B-10 capture, 4 for H-1 capture 
};

/**
 * @brief class of all particles' histories in one nps simulation.
 * 
 */
typedef st_tree::tree<NeutronHistory> NPSHistory;
typedef st_tree::tree<NeutronHistory>::iterator NPSHistoryIterator;
typedef st_tree::tree<NeutronHistory>::const_iterator NPSHistoryConstIterator;