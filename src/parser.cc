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

MCNPPTRACBinary::MCNPPTRACBinary(std::string const &ptracPath) : ptracFile(ptracPath, std::ios_base::binary)
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
    VariableIDNum idnum = VariableIDNum{nbDataNPS, nbDataSrcLong, nbDataSrcDouble,
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

    indices = EventIndices{ 0, // event type
                            5, // cell num
                            0, // x
                            1, // y
                            2, // z
                            6, // energy
                            7, // weight
                            8, // time
                            idnum
                            };
}

void MCNPPTRACBinary::parsePTRACRecord()
{
    npsHistory = NPSHistory(0);
    
    long nps = -1;
    long event = -1, oldEvent = -1;
    constexpr long lastEvent = 9000;
    constexpr long terEvent = 5000;
    constexpr long sourceEvent = 2030;
    long cell = -1;

    double px = 0.;
    double py = 0.;
    double pz = 0.;
    double erg = 0;
    double wt = 0;
    double tme = 0;

    std::string buffer = readRecord(ptracFile); // NPS line
    std::tie(nps, event) = reinterpretBuffer<long, long>(buffer);
    if (event != sourceEvent)
    {
        throw std::logic_error("expected bank event at the start of the history");
    }

    ParticleHistory parHist = ParticleHistory();
    while (event != lastEvent)
    {
        buffer = readRecord(ptracFile); // data line (all doubles, even though the
                                        // first group are actually longs)
        std::stringstream bufferStream(buffer);
        for (int i = 0; i < indices.idNum.nbDataBnkLong; ++i)
        {
            const auto someLong = static_cast<long>(std::get<0>(reinterpretBuffer<double>(bufferStream)));
            if (i == indices.event)
            {
                oldEvent = event;
                event = someLong;
            }
            else if (i == indices.cell)
            {
                cell = someLong;
            }
            // throw away the others
        }
        for (int i = 0; i < indices.idNum.nbDataSrcDouble; ++i)
        {
            const auto someDouble = std::get<0>(reinterpretBuffer<double>(bufferStream));
            if (i == indices.px)
            {
                px = someDouble;
            }
            else if (i == indices.py)
            {
                py = someDouble;
            }
            else if (i == indices.pz)
            {
                pz = someDouble;
            }
            else if (i == indices.erg)
            {
                erg = someDouble;
            }
            else if (i == indices.wt)
            {
                wt = someDouble;
            }
            else if (i == indices.tme)
            {
                tme = someDouble;
            }
            // throw away the others
        }
        parHist.push_back(Event{nps, oldEvent, cell, {px,py,pz}, erg, wt, tme});
        if (event == sourceEvent || event == lastEvent)
        {
            npsHistory.push_back(parHist);
            parHist = ParticleHistory();
        }
    }
}
