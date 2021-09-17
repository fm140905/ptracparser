/**
 * @file pulse.cc
 * @author Ming Fang
 * @brief Extract pulse from particle historie
 * @date 2021-09-13
 */
#include "pulse.hh"
#include "iomanip"

Pulse::Pulse(const ParticleHistory& parHist)
{
    energy = 0;
    pos = std::vector<double>(3, 0);
    for (auto iter = parHist.begin(); iter != parHist.end(); iter++)
    {
        if (iter->cellID == 601 && iter->eventID != 5000)
        {
            if (energy == 0)
            {
                // first event in detector cell
                nps = iter->nps;
                startPos = iter->pos;
                time = iter->time;
            }
            auto nextit = std::next(iter, 1);
            endPos = nextit->pos;
            energy += iter->energy - nextit->energy;
        }
    }
    if (energy > 0)
    {
        for (int i = 0; i < 3; i++)
        {
            pos[i] = (startPos[i] + endPos[i]) * 0.5;
        }
        // getLatticeIndex();
    }
}

/*
void Pulse::getLatticeIndex()
{
    static const double r = 0.25; // radius of straw
    static const double delta = 10.0 / 11.0; //  Distance between two straws, in cm.
    const double x = pos[0];
    const double y = pos[1];
    const double SQRT_3 = 1.73205080757;
    const double i_low = 22 + 2*(y-r) / (SQRT_3 * delta);
    const double i_high = 22 + 2*(y+r) / (SQRT_3 * delta);
    const int i = std::ceil(i_low);
    if (i != int(std::floor(i_high)))
    {
        std::cout << "Cannot find row index. nps = " << nps
                  << ", x = " << x  
                  << ", y = " << y
                  << ",i_low = " << i_low
                  << ", i_high = " << i_high << std::endl;
    }
    double j_low = 22 + (x-r) / delta - (std::ceil(i_low)-22) / 2;
    double j_high = 22 + (x+r) / delta - (std::ceil(i_low)-22) / 2;
    const int j = std::ceil(j_low);
    if (j != int(std::floor(j_high)))
    {
        std::cout << "Cannot find column index. nps = " << nps
                  << ", x = " << x 
                  << ", y = " << y
                  << ",j_low = " << j_low
                  << ", j_high = " << j_high << std::endl;
    }
    latticeIndex = std::make_pair(i, j);
}
*/
std::ostream &operator<<(std::ostream &os, const Pulse& p)
{
    os << std::fixed << std::setprecision(6)
       << std::setw(12) << p.startPos[0]
       << std::setw(12) << p.startPos[1]
       << std::setw(12) << p.startPos[2]
       << std::setw(12) << p.endPos[0]
       << std::setw(12) << p.endPos[1]
       << std::setw(12) << p.endPos[2]
       << std::setw(12) << p.energy
       << std::setw(24) << p.time
    //    << std::setw(8) << p.latticeIndex.first
    //    << std::setw(8) << p.latticeIndex.second
       << std::setw(12) << p.nps;
    return os;
}