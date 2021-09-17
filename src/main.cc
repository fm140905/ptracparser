/*
 * @Description: Main file for processing ptrac file.
 * @Author: Ming Fang
 * @Date: 2021-09-14 11:43:25
 * @LastEditors: Ming Fang
 * @LastEditTime: 2021-09-14 12:54:21
 */
#include <fstream>
#include <string>
#include <chrono>

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

    const std::string ptracFilePath(argv[1]);
    MCNPPTRACBinary ptracFile(ptracFilePath);
    const int maxNum(1e6);
    int pulseIdx(0);
    std::vector<Pulse> pulses;
    // read pulses from file
    int nps(0);
    while (ptracFile.readNextNPS(1e6))
    {
        if (ptracFile.getNPSRead() % 1000 == 0)
        {
            std::cout << "NPS = " << ptracFile.getNPSRead() << '\n';
        }
        const NPSHistory record = ptracFile.getNPSHistory();
        for (auto iter = record.begin(); iter != record.end(); iter++)
        {
            Pulse newPulse(*iter);
            if (newPulse.energy <= 0)
                continue;
            // else, it is a vaild pulse
            pulses.push_back(newPulse);
            pulseIdx++;
            if (pulseIdx >= maxNum)
                break;
        }
        if (pulseIdx >= maxNum)
            break;
    }

    // write header
    outfile << "#    x1(cm)      y1(cm)      z1(cm)      x2(cm)      y2(cm)      z2(cm)    energy(MeV)        time(shakes)           nps\n";
    // write pulses to file
    for (int i = 0; i < pulses.size(); i++)
    {
        outfile << pulses[i] << '\n';
    }
    outfile.close();

    auto endTime = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << "ms" << std::endl; 
    
    return 0;
}