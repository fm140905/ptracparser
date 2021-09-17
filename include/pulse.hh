/**
 * @file pulse.hh
 * @author Ming Fang
 * @brief Generate pulse based on particle history
 * @date 2021-09-13
 */
#pragma once
#include "parser.hh"
#include <math.h>

class Pulse
{
    friend std::ostream &operator<<(std::ostream &os, const Pulse& p);
private:
    /**
     * @brief Get the Lattice Index of the cell where the pulse is created
     * 
     */
    // void getLatticeIndex();
public:
    Pulse() = default;
    Pulse(const ParticleHistory& parHist);
    /* data */
    long nps; // particle nps
    // std::pair<int, int> latticeIndex;
    std::vector<double> startPos; // start of track, (x,y,z), cm, for debugging
    std::vector<double> endPos; // end of track, (x,y,z), cm, for debugging
    std::vector<double> pos; // center of track, (x,y,z), cm
    double time; // time stamp, shakes
    double energy; // deposited energy, MeV
};

std::ostream &operator<<(std::ostream &os, const Pulse& p);



