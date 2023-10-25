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
#include <iomanip>
#include <algorithm>

#include "parser.hh"
#include "pulse.hh"
int main(int argc, char** argv)
{
    auto startTime = std::chrono::high_resolution_clock::now();

    const std::string outpath("pulses.txt");
    std::ofstream outfile;
    outfile.open(outpath, std::ios::out);
    if (!outfile.good())
    {
        throw std::invalid_argument("Cannot create file: " + outpath);
    }

    // const std::string ptracFilePath(argv[1]);
    const std::string ptracFilePath("/media/ming/Elements/multiplicity_test/BCS_NMC_1500/ptrac_parser/ptracparser/tests/testdata/v6/ptrac");
    MCNPPTRACBinary ptracFile(ptracFilePath);
    const long maxNum(1e9);
    long pulseIdx(0);
    std::vector<double> pulse_timestamps;
    // read pulses from file
    // int nps(0);
    while (ptracFile.readNextNPS(maxNum))
    {
        if (ptracFile.getNPSRead() % 1000000 == 0)
        {
            std::cout << "NPS = " << ptracFile.getNPSRead() << '\n';
        }
        const NPSHistory record = ptracFile.getNPSHistory();
        for (auto iter = record.begin(); iter != record.end(); iter++)
        {
            for (auto iter2 = iter->cbegin(); iter2 != iter->cend(); iter2++)
            {
                // if cellID != 100, continue
                if (iter2->cellID != 100)
                    continue;
                // else, it is a valid event
                pulse_timestamps.push_back(iter2->time);
                pulseIdx++;
                if (pulseIdx >= maxNum)
                    break;
            }
            if (pulseIdx >= maxNum)
                break;
        }
        if (pulseIdx >= maxNum)
            break;
    }

    // sort timestamps
    std::sort(pulse_timestamps.begin(), pulse_timestamps.end());
    // write header
    outfile << "#    time(shakes)\n";
    // write pulses to file
    for (int i = 0; i < pulse_timestamps.size(); i++)
    {
        outfile << std::fixed << std::setprecision(8) << pulse_timestamps[i] << '\n';
    }
    outfile.close();

    auto endTime = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << "ms" << std::endl; 
    
    return 0;
}