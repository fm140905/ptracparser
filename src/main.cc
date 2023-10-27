/*
 * @Description: Main file for processing ptrac file.
 * @Author: Ming Fang
 * @Date: 2021-09-14 11:43:25
 * @LastEditors: Ming Fang
 * @LastEditTime: 2022-07-28 00:07:54
 */
#include <fstream>
#include <string>
#include <chrono>

// #include <highfive/H5DataSet.hpp>
// #include <highfive/H5DataSpace.hpp>
// #include <highfive/H5File.hpp>

#include "parser.hh"
int main(int argc, char** argv)
{
    auto startTime = std::chrono::high_resolution_clock::now();

    // const std::string ptracFilePath(argv[1]);
    const std::string ptracFilePath("/media/ming/Elements/multiplicity_test/BCS_NMC_1500/cf252_sphere/ptrac/v14/ptrac");
    MCNPPTRACBinary ptracFile(ptracFilePath);
    const long maxNum(10000);
    // std::vector<ParticleHistory> spontaneous_fission_neutrons;
    // std::vector<ParticleHistory> induced_fission_neutrons;

    bool flag = true;
    std::list<NPSHistory> nps_history;
    while (ptracFile.readNextNPS(maxNum))
    {
        
        if (ptracFile.getNPSRead() % 10000 == 0)
        {
            std::cout << "NPS = " << ptracFile.getNPSRead() << '\n';
        }
        const NPSHistory record = ptracFile.getNPSHistory();
        nps_history.push_back(record);
    }

    // extract SF multiplicity distribution
    // extract IF multiplicity distribution
    // extract number of SF reactions, number of IF reactions, number of capture reactions
    // extract number of neutrons created in SF, number of neutrons created in IF, number of neutrons terminated in capture
    // extract number of neutrons terminated in B-10 capture
    // extract the time of B-10 capture
    // extract the time difference between the creation and termination of neutrons

    auto endTime = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << "ms" << std::endl; 
    
    return 0;
}