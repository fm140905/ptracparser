/**
* @file Parser.cc
*
*
* @brief Binary PTRAC file parser
*
* @author J. Faustin, Ming Fang
* @version 1.1
*/

#include "parser.hh"
#include <algorithm>
#include <cassert>
#include <cctype>
#include <unistd.h>

/************************************
*                                  *
*  methods of the MCNPPTRAC class  *
*                                  *
************************************/

MCNPPTRAC::MCNPPTRAC() : npsRead(0)
{
}

void MCNPPTRAC::incrementNPSRead()
{
    ++npsRead;
}

long MCNPPTRAC::getNPSRead()
{
    return npsRead;
}

NPSHistory const &MCNPPTRAC::getNPSHistory() const
{
    return npsHistory;
}

/*****************************************
*                                        *
*  methods of the MCNPPTRACBinary class  *
*                                        *
******************************************/

MCNPPTRACBinary::MCNPPTRACBinary(std::string const &ptracPath)
 : ptracFile(ptracPath, std::ios_base::binary)
{
    if (ptracFile.fail())
    {
        std::cerr << "PTRAC file " << ptracPath << " not found." << std::endl;
        exit(EXIT_FAILURE);
    }
    parseHeader();
}

bool MCNPPTRACBinary::readNextNPS(long maxReadHist)
{
    if ((ptracFile && ptracFile.peek() != EOF) && npsRead <= maxReadHist)
    {
        parsePTRACRecord();
        incrementNPSRead();
        return true;
    }
    return false;
}

void MCNPPTRACBinary::parseHeader()
{
    skipHeader();
    skipPtracInputData();
    parseVariableIDs();
}

void MCNPPTRACBinary::skipHeader()
{
    readRecord(ptracFile); // header
    readRecord(ptracFile); // code, version, dates
    readRecord(ptracFile); // calculation title
}

void MCNPPTRACBinary::skipPtracInputData()
{
    std::string buffer = readRecord(ptracFile);
    std::stringstream bufferStream(buffer);
    int n_fields_total = (int)readBinary<double>(bufferStream);
    int n_fields_read = 0;
    int n_desc = 0;

    while (n_fields_read < n_fields_total)
    {
        if (bufferStream.peek() == EOF)
        {
            buffer = readRecord(ptracFile);
            bufferStream.str(buffer);
            bufferStream.clear();
        }
        if (n_desc == 0)
        {
            n_desc = (int)readBinary<double>(bufferStream);
            ++n_fields_read;
        }
        else
        {
            --n_desc;
            readBinary<double>(bufferStream); // throw away field descriptor
        }
    }
}

void MCNPPTRACBinary::parseVariableIDs()
{
    std::string buffer = readRecord(ptracFile); // line 6
    std::stringstream bufferStream(buffer);
    //   auto const fields = reinterpretBuffer<int, long, long>(buffer);
    // number of variables on NPS line
    const int nbDataNPS = readBinary<int>(bufferStream);
    // number of variable on first line for an src event
    const long nbDataSrcLong = readBinary<long>(bufferStream);
    // number of variable on second line for an src event
    const long nbDataSrcDouble = readBinary<long>(bufferStream);
    // number of variable on first line for an bank event
    const long nbDataBnkLong = readBinary<long>(bufferStream);
    // number of variable on second line for an bank event
    const long nbDataBnkDouble = readBinary<long>(bufferStream);
    // number of variable on first line for an suf event
    const long nbDataSufLong = readBinary<long>(bufferStream);
    // number of variable on second line for an suf event
    const long nbDataSufDouble = readBinary<long>(bufferStream);
    // number of variable on first line for an col event
    const long nbDataColLong = readBinary<long>(bufferStream);
    // number of variable on second line for an col event
    const long nbDataColDouble = readBinary<long>(bufferStream);
    // number of variable on first line for an ter event
    const long nbDataTerLong = readBinary<long>(bufferStream);
    // number of variable on second line for an ter event
    const long nbDataTerDouble = readBinary<long>(bufferStream);
    idNum = VariableIDNum{nbDataNPS, nbDataSrcLong, nbDataSrcDouble,
                                        nbDataBnkLong, nbDataBnkDouble,
                                        nbDataSufLong, nbDataSufDouble,
                                        nbDataColLong, nbDataColDouble,
                                        nbDataTerLong, nbDataTerDouble};

    buffer = readRecord(ptracFile); // line 7
    bufferStream.str(buffer);
    bufferStream.clear();

    // Skip over the NPS data line. For some reason these fields are written as
    // longs.
    for (int i = 0; i < nbDataNPS; ++i)
    {
        reinterpretBuffer<long>(bufferStream); // variable ids on NPS lines
    }

    // throw away variable id on event lines, int
}

void MCNPPTRACBinary::parsePTRACRecord()
{
    npsHistory = NPSHistory(0);
    
    long nextEvent = -1, currEvent = -1;
    constexpr long lastEvent = 9000;

    std::vector<long> firstGroup(idNum.nbDataBnkLong);
    std::vector<double> secondGroup(idNum.nbDataBnkDouble);

    long nps;
    std::string buffer = readRecord(ptracFile); // NPS line
    std::tie(nps, nextEvent) = reinterpretBuffer<long, long>(buffer);
    if (nextEvent != 1000)
    {
        std::cout << "Event number: " << nextEvent << std::endl;
        throw std::logic_error("expected source event at the start of the history");
    }

    ParticleHistory parHist{Event(), Event()};
    while (nextEvent != lastEvent)
    {
        currEvent = nextEvent;
        buffer = readRecord(ptracFile); // data line (all doubles, even though the
                                        // first group are actually longs)
        std::stringstream bufferStream(buffer);
        switch (currEvent)
        {
        case 1000:
            // spontaneous fission neutron creation
            parseBuffer(bufferStream, firstGroup, secondGroup, idNum.nbDataSrcLong);
            // select fields that we need
            nextEvent = firstGroup[indices.event];
            parHist[0].nps = nps;
            parHist[0].eventID = currEvent;
            parHist[0].pos = {secondGroup[indices.px],
                                   secondGroup[indices.py],
                                   secondGroup[indices.pz]};
            parHist[0].energy = secondGroup[indices.erg];
            parHist[0].time = secondGroup[indices.tme];
            // debugging
// #ifndef NDEBUG
            // if third element of first group != 40, then throw an error
            if (firstGroup[2] != 40)
            {
                std::string message = "NPS="+std::to_string(nps) + ". Expected source type 40, but got " + std::to_string(firstGroup[2]);
                std::cout << message << std::endl; 
                // throw std::logic_error("NPS="+std::to_string(nps) + ". Expected source type 40, but got " + std::to_string(firstGroup[2]));
            }
            // if the cell number (fourth element) is not 100, then throw an error
            if (firstGroup[3] != 100)
            {
                std::string message = "NPS="+std::to_string(nps) + ". Expected cell number 100, but got " + std::to_string(firstGroup[3]);
                std::cout << message << std::endl;
                // throw std::logic_error("NPS="+std::to_string(nps) + ". Expected cell number 100, but got " + std::to_string(firstGroup[indices.cell]));
            }
// #endif
            break;
        case 2000:
            // spontaneous fission neutron creation
            parseBuffer(bufferStream, firstGroup, secondGroup, idNum.nbDataBnkLong);
            // select fields that we need
            nextEvent = firstGroup[indices.event];
            parHist[0].nps = nps;
            parHist[0].eventID = currEvent;
            parHist[0].pos = {secondGroup[indices.px],
                                   secondGroup[indices.py],
                                   secondGroup[indices.pz]};
            parHist[0].energy = secondGroup[indices.erg];
            parHist[0].time = secondGroup[indices.tme];

            // debugging
// #ifndef NDEBUG
            // if the third element of first group != 98252, then throw an error
            if (firstGroup[2] != 98252)
            {
                std::string message = "NPS="+std::to_string(nps) + ". Expected nuclide 98252, but got " + std::to_string(firstGroup[2]);
                std::cout << message << std::endl;
                // throw std::logic_error("NPS="+std::to_string(nps) + ". Expected nuclide 98252, but got " + std::to_string(firstGroup[2]));
            }
            // if the fourth element of first group != 0, then throw an error
            if (firstGroup[3] != 0)
            {
                std::string message = "NPS="+std::to_string(nps) + ". Expected reaction type 0, but got " + std::to_string(firstGroup[3]);
                std::cout << message << std::endl;
                // throw std::logic_error("NPS="+std::to_string(nps) + ". Expected reaction type 0, but got " + std::to_string(firstGroup[3]));
            }
            // if the cell number is not 100, then throw an error
            if (firstGroup[indices.cell] != 100)
            {
                std::string message = "NPS="+std::to_string(nps) + ". Expected cell number 100, but got " + std::to_string(firstGroup[indices.cell]);
                std::cout << message << std::endl;
                // throw std::logic_error("NPS="+std::to_string(nps) + ". Expected cell number 100, but got " + std::to_string(firstGroup[indices.cell]));
            }
// #endif
            break;
        case 2007:
            // induced fission neutron creation
            parseBuffer(bufferStream, firstGroup, secondGroup, idNum.nbDataBnkLong);
            // select fields that we need
            nextEvent = firstGroup[indices.event];
            parHist[0].nps = nps;
            parHist[0].eventID = currEvent;
            parHist[0].pos = {secondGroup[indices.px],
                                   secondGroup[indices.py],
                                   secondGroup[indices.pz]};
            parHist[0].energy = secondGroup[indices.erg];
            parHist[0].time = secondGroup[indices.tme];

            // debugging
// #ifndef NDEBUG
            // if the third element of first group != 98252, then throw an error
            if (firstGroup[2] != 98252)
            {
                std::string message = "NPS="+std::to_string(nps) + ". Expected nuclide 98252, but got " + std::to_string(firstGroup[2]);
                std::cout << message << std::endl;
                // throw std::logic_error("NPS="+std::to_string(nps) + ". Expected nuclide 98252, but got " + std::to_string(firstGroup[2]));
            }
            // if the fourth element of first group != 19, then throw an error
            if (firstGroup[3] != 19)
            {
                std::string message = "NPS="+std::to_string(nps) + ". Expected reaction type 19, but got " + std::to_string(firstGroup[3]);
                std::cout << message << std::endl;
                // throw std::logic_error("NPS="+std::to_string(nps) + ". Expected reaction type 19, but got " + std::to_string(firstGroup[3]));
            }
            // if the cell number is not 100, then throw an error
            if (firstGroup[indices.cell] != 100)
            {
                std::string message = "NPS="+std::to_string(nps) + ". Expected cell number 100, but got " + std::to_string(firstGroup[indices.cell]);
                std::cout << message << std::endl;
                // throw std::logic_error("NPS="+std::to_string(nps) + ". Expected cell number 100, but got " + std::to_string(firstGroup[indices.cell]));
            }
// #endif
            break;
        case 5000:
            // termination
            parseBuffer(bufferStream, firstGroup, secondGroup, idNum.nbDataTerLong);
            // select fields that we need
            nextEvent = firstGroup[indices.event];
            parHist[1].nps = nps;
            parHist[1].eventID = currEvent;
            parHist[1].pos = {secondGroup[indices.px],
                                   secondGroup[indices.py],
                                   secondGroup[indices.pz]};
            parHist[1].energy = secondGroup[indices.erg];
            parHist[1].time = secondGroup[indices.tme];

            // debugging
// #ifndef NDEBUG
            // if the third element of first group != 12, then throw an error
            if (firstGroup[2] != 12)
            {
                std::string message = "NPS="+std::to_string(nps) + ". Expected termination type 12, but got " + std::to_string(firstGroup[2]);
                std::cout << message << std::endl;
                // throw std::logic_error("NPS="+std::to_string(nps) + ". Expected termination type 12, but got " + std::to_string(firstGroup[2]));
            }
            // if the cell number is not 100, then throw an error
            if (firstGroup[indices.cell] != 100)
            {
                std::string message = "NPS="+std::to_string(nps) + ". Expected cell number 100, but got " + std::to_string(firstGroup[indices.cell]);
                std::cout << message << std::endl;
                // throw std::logic_error("NPS="+std::to_string(nps) + ". Expected cell number 100, but got " + std::to_string(firstGroup[indices.cell]));
            }
// #endif
            // add the parHist to npsHistory
            npsHistory.push_back(parHist);
           break;
        default:
            parseBuffer(bufferStream, firstGroup, secondGroup, idNum.nbDataBnkLong);
            std::string message = "NPS="+std::to_string(nps) + ". Unexpected event type " + std::to_string(currEvent);
            std::cout << message << std::endl;
            // throw std::logic_error("NPS="+std::to_string(nps) + ". Unexpected event type" + std::to_string(currEvent));
            break;
        }
    }
}

void MCNPPTRACBinary::parseBuffer(std::stringstream& bufferStream, std::vector<long>& firstGroup,
    std::vector<double>& secondGroup, int size)
{
    for (int i = 0; i < size; i++)
    {
        const auto someLong = static_cast<long>(std::get<0>(reinterpretBuffer<double>(bufferStream)));
        firstGroup[i] = someLong;
    }
    for (int i = 0; i < secondGroup.size(); i++)
    {
        const auto someDouble = std::get<0>(reinterpretBuffer<double>(bufferStream));
        secondGroup[i] = someDouble;
    }
}