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

MCNPPTRAC::MCNPPTRAC() : nbPointsRead(0)
{
}

void MCNPPTRAC::incrementNbPointsRead()
{
    ++nbPointsRead;
}

long MCNPPTRAC::getNbPointsRead()
{
    return nbPointsRead;
}

History const &MCNPPTRAC::getPTRACRecord() const
{
    return history;
}

// /*****************************************
// *                                       *
// *  methods of the MCNPPTRACASCII class  *
// *                                       *
// *****************************************/

// MCNPPTRACASCII::MCNPPTRACASCII(std::string const &ptracPath) : MCNPPTRAC(),
//                                                                nbDataCellMaterialLine(0),
//                                                                ptracFile(ptracPath)
// {
//     if (ptracFile.fail())
//     {
//         std::cerr << "PTRAC file " << ptracPath << " not found." << std::endl;
//         exit(EXIT_FAILURE);
//     }
//     // The number of header lines must be 8!!
//     goThroughHeaderPTRAC(8);
// }

// void MCNPPTRACASCII::getNextLinePtrac()
// {
//     getline(ptracFile, currentLine);
// }

// std::pair<int, int> MCNPPTRACASCII::readPointEvent()
// {
//     std::istringstream iss(currentLine);
//     int pointID, eventID;
//     iss >> pointID >> eventID;
//     return {pointID, eventID};
// }

// std::pair<int, int> MCNPPTRACASCII::readCellMaterial()
// {
//     std::istringstream iss(currentLine);
//     int dummy, volumeID, materialID;
//     for (int ii = 1; ii <= nbDataCellMaterialLine - 2; ii++)
//     {
//         iss >> dummy;
//     }
//     iss >> volumeID >> materialID;
//     return {volumeID, materialID};
// }

// std::vector<double> MCNPPTRACASCII::readPoint()
// {
//     std::istringstream iss(currentLine);
//     double pointX, pointY, pointZ;
//     iss >> pointX >> pointY >> pointZ;
//     return {pointX, pointY, pointZ};
// }

// bool MCNPPTRACASCII::readNextPtracData(long maxReadPoint)
// {
//     if ((ptracFile && !ptracFile.eof()) && (getNbPointsRead() <= maxReadPoint))
//     {
//         getline(ptracFile, currentLine);
//         if (!currentLine.empty())
//         {
//             auto const pointEvent = readPointEvent();
//             getline(ptracFile, currentLine);
//             auto const cellMaterial = readCellMaterial();
//             getline(ptracFile, currentLine);
//             auto const point = readPoint();
//             incrementNbPointsRead();
//             // TODO
//             //   record = {pointEvent.first, pointEvent.second,
//             //             cellMaterial.first, cellMaterial.second,
//             //             point};
//             return true;
//         }
//         else
//         {
//             return false;
//         }
//     }
//     else
//     {
//         return false;
//     }
// }

// void MCNPPTRACASCII::goThroughHeaderPTRAC(int nHeaderLines)
// {
//     std::string line5, line6;
//     for (int ii = 0; ii < nHeaderLines; ii++)
//     {
//         getline(ptracFile, currentLine);
//         if (ii == 5)
//         {
//             line5 = currentLine;
//         }
//         if (ii == 6)
//         {
//             line6 = currentLine;
//         }
//     }
//     int nbData = getDataFromLine5Ptrac(line5);
//     checkDataFromLine6Ptrac(line6, nbData);
// }

// int MCNPPTRACASCII::getDataFromLine5Ptrac(const std::string &line5)
// {
//     int nbDataPointEventLine;
//     std::istringstream iss(line5);
//     iss >> nbDataPointEventLine >> nbDataCellMaterialLine;
//     int nbData = nbDataPointEventLine + nbDataCellMaterialLine;
//     return nbData;
// }

// void MCNPPTRACASCII::checkDataFromLine6Ptrac(const std::string &line6, int nbData)
// {
//     std::vector<int> data(nbData);
//     constexpr int cellIDPtracCode = 17;
//     constexpr int materialIDPtracCode = 18;
//     std::istringstream iss(line6);
//     for (int jj = 0; jj < nbData; jj++)
//     {
//         iss >> data[jj];
//     }
//     if ((data[nbData - 2] != cellIDPtracCode) && (data[nbData - 1] != materialIDPtracCode))
//     {
//         std::cerr << "PTRAC file format not suitable. Please see Oracle/data/slapb file for example..." << std::endl;
//         exit(EXIT_FAILURE);
//     }
// }

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

bool MCNPPTRACBinary::readNextPtracData(long maxReadPoint)
{
    if ((ptracFile && ptracFile.peek() != EOF) && getNbPointsRead() <= maxReadPoint)
    {
        parsePTRACRecord();
        incrementNbPointsRead();
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

    indices = PTRACRecordIndices { 0, // event type
                                   5, // cell num
                                   0, // x
                                   1, // y
                                   2, // z
                                   6, // energy
                                   7, // weight
                                   8 // time
                                   };
    //   constexpr int eventID = 7;
    //   constexpr int cellID = 17;
    //   constexpr int
    //   int eventIdx = -1;
    //   int cellIdx = -1;
    //   int matIdx = -1;
    //   for (int i = 0; i < nbDataSrcLong; ++i) {
    //     const int varID = std::get<0>(reinterpretBuffer<int>(bufferStream));
    //     switch (varID) {
    //     case eventID:
    //       eventIdx = i;
    //       break;
    //     case cellID:
    //       cellIdx = i;
    //       break;
    //     case matID:
    //       matIdx = i;
    //       break;
    //     }
    //   }

    //   constexpr int pointXID = 20;
    //   constexpr int pointYID = 21;
    //   constexpr int pointZID = 22;
    //   int pointXIdx = -1;
    //   int pointYIdx = -1;
    //   int pointZIdx = -1;
    //   for (int i = 0; i < nbDataSrcDouble; ++i) {
    //     const int varID = std::get<0>(reinterpretBuffer<int>(bufferStream));
    //     switch (varID) {
    //     case pointXID:
    //       pointXIdx = i;
    //       break;
    //     case pointYID:
    //       pointYIdx = i;
    //       break;
    //     case pointZID:
    //       pointZIdx = i;
    //       break;
    //     }
    //   }

    //   indices = PTRACRecordIndices{eventIdx, cellIdx, matIdx, pointXIdx, pointYIdx, pointZIdx, nbDataSrcLong, nbDataSrcDouble};
}

void MCNPPTRACBinary::parsePTRACRecord()
{
    history = History();
    
    long point = -1;
    long event = -1, oldEvent = -1;
    constexpr long lastEvent = 9000;
    constexpr long sourceEvent = 2030;
    long cell = -1;

    double px = 0.;
    double py = 0.;
    double pz = 0.;
    double erg = 0;
    double wt = 0;
    double tme = 0;

    std::string buffer = readRecord(ptracFile); // NPS line
    std::tie(point, event) = reinterpretBuffer<long, long>(buffer);
    if (event != sourceEvent)
    {
        throw std::logic_error("expected bank event at the start of the history");
    }

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
        history.push_back(PTRACRecord{point, oldEvent, cell, {px,py,pz}, erg, wt, tme});
    }
}
