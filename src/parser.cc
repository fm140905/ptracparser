/**
* @file Parser.cc
*
*
* @brief Binary PTRAC file parser
*
* @author J. Faustin, Ming Fang
* @version 1.1
*/

#include <algorithm>
#include <cassert>
#include <cctype>
#include <unistd.h>
#include <iomanip>

#include "parser.hh"
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
    npsHistory = NPSHistory();
    // create a virtual root event representing the simultaneous creation of SF neutrons
    NeutronHistory root_event = NeutronHistory();
    
    long nextEvent = -1, currEvent = -1;
    int next_creation_reaction = -1;
    int next_destruction_reaction = -1;
    constexpr long lastEvent = 9000;

    std::vector<long> firstGroup(idNum.nbDataBnkLong);
    std::vector<double> secondGroup(idNum.nbDataBnkDouble);

    long nps;
    std::string buffer = readRecord(ptracFile); // NPS line
    std::tie(nps, nextEvent) = reinterpretBuffer<long, long>(buffer);
    root_event.nps = nps;
    if (nextEvent != 1000)
    {
        std::cout << "Event number: " << nextEvent << std::endl;
        throw std::logic_error("expected source event at the start of the history");
    }
// #ifndef NDEBUG
//     // print nps and nextEvent
//     if (nps == 145234)
//         std::cout << "NPS: " << nps << std::endl;
// #endif

    NeutronHistory curr_neutron = NeutronHistory();
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
            parseBuffer(bufferStream, firstGroup, secondGroup, idNum.nbDataSrcLong, nps);
            // update root event
            root_event.creation_reaction = -1; 
            root_event.creation_pos = {secondGroup[indices.px],
                                   secondGroup[indices.py],
                                   secondGroup[indices.pz]};
            root_event.creation_energy = -1;
            root_event.creation_time = secondGroup[indices.tme];
            root_event.destruction_reaction = SF_REACTION;
            root_event.destruction_pos = root_event.creation_pos;
            root_event.destruction_energy = -1;
            root_event.destruction_time = root_event.creation_time;
            npsHistory.insert(root_event);

            // select fields that we need
            nextEvent = firstGroup[indices.event];
            curr_neutron.nps = nps;
            curr_neutron.creation_reaction = SF_REACTION;
            curr_neutron.creation_pos = root_event.creation_pos;
            curr_neutron.creation_energy = secondGroup[indices.erg];
            curr_neutron.creation_time = secondGroup[indices.tme];
            // debugging
#ifndef NDEBUG
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
#endif
            break;
        case 2000:
            // spontaneous fission neutron creation
            parseBuffer(bufferStream, firstGroup, secondGroup, idNum.nbDataBnkLong, nps);
            // select fields that we need
            nextEvent = firstGroup[indices.event];
            curr_neutron.creation_reaction = SF_REACTION;
            curr_neutron.creation_pos = {secondGroup[indices.px],
                                   secondGroup[indices.py],
                                   secondGroup[indices.pz]};
            curr_neutron.creation_energy = secondGroup[indices.erg];
            curr_neutron.creation_time = secondGroup[indices.tme];

            // debugging
#ifndef NDEBUG
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
#endif
            break;
        case 2007:
            // induced fission neutron creation or (n,2n) neutron creation
            parseBuffer(bufferStream, firstGroup, secondGroup, idNum.nbDataBnkLong, nps);
            // select fields that we need
            nextEvent = firstGroup[indices.event];
            curr_neutron.creation_reaction = next_creation_reaction;
            curr_neutron.creation_pos = {secondGroup[indices.px],
                                   secondGroup[indices.py],
                                   secondGroup[indices.pz]};
            curr_neutron.creation_energy = secondGroup[indices.erg];
            curr_neutron.creation_time = secondGroup[indices.tme];

            // debugging
#ifndef NDEBUG
            // if the third element of first group != 98252, then throw an error
            if (firstGroup[2] != 98252)
            {
                std::string message = "NPS="+std::to_string(nps) + ". Expected nuclide 98252, but got " + std::to_string(firstGroup[2]);
                std::cout << message << std::endl;
                // throw std::logic_error("NPS="+std::to_string(nps) + ". Expected nuclide 98252, but got " + std::to_string(firstGroup[2]));
            }
            // if the fourth element of first group != 19 or 2 or 3, then throw an error
            if (firstGroup[3] != 19 && firstGroup[3] != 2 && firstGroup[3] != 3 && firstGroup[3] != 4)
            {
                std::string message = "NPS="+std::to_string(nps) + ". Expected reaction type 19 or 2 or 3 or 4, but got " + std::to_string(firstGroup[3]);
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
#endif
            break;
        case 4000:
            // collision
            parseBuffer(bufferStream, firstGroup, secondGroup, idNum.nbDataColLong, nps);
            // select fields that we need
            nextEvent = firstGroup[indices.event];
            // if scattering with 1001/5010, do nothing
            static const std::map<std::pair<long, long>, int> special_collisions = {{{98252, 18}, IF_REACTION}, // fission
                                                                                    {{98252, 16}, N2N_REACTION}, // (n,2n)
                                                                                    {{98252, 17}, N3N_REACTION}, // (n, 3n)
                                                                                    {{98252, 37}, N4N_REACTION}}; // (n, 4n)
            // if pair is in special_collisions, then set next_creation_reaction to the value
            if (special_collisions.count({firstGroup[2], firstGroup[3]}) > 0)
            {
                next_creation_reaction = special_collisions.at({firstGroup[2], firstGroup[3]});
                curr_neutron.destruction_reaction = next_creation_reaction;
                curr_neutron.destruction_pos = {secondGroup[indices.px],
                                    secondGroup[indices.py],
                                    secondGroup[indices.pz]};
                curr_neutron.destruction_energy = 0;
                curr_neutron.destruction_time = secondGroup[indices.tme];
                // find parent neutron based on pos and tme
                // add curr_neutron to npsHistory
                _add_neutron(curr_neutron);

                curr_neutron.creation_reaction = next_creation_reaction;
                curr_neutron.creation_pos = curr_neutron.destruction_pos;
                curr_neutron.creation_energy = secondGroup[indices.erg];
                curr_neutron.creation_time = curr_neutron.destruction_time;
            }
            else
            {
                if ( !(firstGroup[3] == 2) && !(firstGroup[2]==5010 && firstGroup[3] >= 50 && firstGroup[3] <= 91)
                                            && !(firstGroup[2]==98252 && firstGroup[3] == 91)) // 
                {
                    std::string message = "NPS="+std::to_string(nps) + ". Unexpected collision type " + std::to_string(firstGroup[2]) + " " + std::to_string(firstGroup[3]);
                    std::cout << message << std::endl;
                }
            }

            // next event is termination
            if (nextEvent == 5000)
            {
                switch (firstGroup[2])
                {
                case 5010:
                    // B-10 capture
                    next_destruction_reaction = B10_CAPTURE;
                    break;
                case 1001:
                    // H-1 capture
                    next_destruction_reaction = H1_CAPTURE;
                    break;
                case 98252:
                    // CF-252 capture
                    next_destruction_reaction = CF252_CAPTURE;
                    break;
                default:
                    // unexpected destruction type.
                    std::string message = "NPS="+std::to_string(nps) + ". Unexpected destruction type " + std::to_string(firstGroup[2]);
                    std::cout << message << std::endl;
                    break;
                }
            }
            break;
        case 5000:
            // termination
            parseBuffer(bufferStream, firstGroup, secondGroup, idNum.nbDataTerLong, nps);
            // select fields that we need
            nextEvent = firstGroup[indices.event];
            curr_neutron.destruction_reaction = next_destruction_reaction;
            if (firstGroup[2] == 14)
            {
                curr_neutron.destruction_reaction = IF_REACTION; // IF with zero neutron emission
            }
            curr_neutron.destruction_pos = {secondGroup[indices.px],
                                   secondGroup[indices.py],
                                   secondGroup[indices.pz]};
            curr_neutron.destruction_energy = secondGroup[indices.erg];
            curr_neutron.destruction_time = secondGroup[indices.tme];
            // debugging
#ifndef NDEBUG
            // if the third element of first group != 12, then throw an error
            if (firstGroup[2] != 12 && firstGroup[2] != 14)
            {
                std::string message = "NPS="+std::to_string(nps) + ". Expected termination type 12 or 14, but got " + std::to_string(firstGroup[2]);
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
#endif
            // add curr_neutron to npsHistory
            _add_neutron(curr_neutron);
           break;
        default:
            parseBuffer(bufferStream, firstGroup, secondGroup, idNum.nbDataBnkLong, nps);
            std::string message = "NPS="+std::to_string(nps) + ". Unexpected event type " + std::to_string(currEvent);
            std::cout << message << std::endl;
            // throw std::logic_error("NPS="+std::to_string(nps) + ". Unexpected event type" + std::to_string(currEvent));
            break;
        }
    }
}

void MCNPPTRACBinary::parseBuffer(std::stringstream& bufferStream, std::vector<long>& firstGroup,
    std::vector<double>& secondGroup, int size, long nps)
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
#ifndef NDEBUG
    // print first group in the first line, second group in the second line. use six digits after decimal point for double
    if (nps == 145234)
    {
        for (auto i : firstGroup)
        {
            std::cout << std::setw(10) << i << " ";
        }
        std::cout << std::endl;
        for (auto i : secondGroup)
        {
            std::cout << std::setw(10) << std::setprecision(6) << std::fixed << i << " ";
        }
        std::cout << std::endl;
    }
#endif
}

void MCNPPTRACBinary::_add_neutron(const NeutronHistory &neutron)
{
    if (neutron.creation_reaction == SF_REACTION)
    {
        // insert to root
        npsHistory.root().insert(neutron);
    }
    else
    {
        // find parent neutron
        for (NPSHistoryIterator j(npsHistory.begin()); j != npsHistory.end(); ++j)
        {
            if (j->data().destruction_time == neutron.creation_time && 
                j->data().destruction_pos[0] == neutron.creation_pos[0] &&
                j->data().destruction_pos[1] == neutron.creation_pos[1] &&
                j->data().destruction_pos[2] == neutron.creation_pos[2])
            {
                j->insert(neutron);
                return;
            }
        }
        // throw error is parent neutron not found
        std::string message = "NPS="+std::to_string(npsHistory.root().data().nps) + ". Parent neutron not found";
        throw std::logic_error(message);
    }
}