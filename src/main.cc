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
#include <numeric>

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

#include "parser.hh"
#include "fission_chain_data_retrieval.hh"

int main(int argc, char** argv)
{
    auto startTime = std::chrono::high_resolution_clock::now();

    // const std::string ptracFilePath(argv[1]);
    const std::string ptracFilePath("/media/ming/Elements/multiplicity_test/BCS_NMC_1500/cf252_sphere/ptrac/v15/ptrac");
    MCNPPTRACBinary ptracFile(ptracFilePath);
    const long maxNum(200000);
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

    // SF multiplicity distribution
    std::vector<int> SF_multiplicities;
    // extract IF multiplicity distribution
    std::vector<int> IF_multiplicities;
    // extract N2N multiplicity distribution
    std::vector<int> N2N_multiplicities;
    // extract N3N multiplicity distribution
    std::vector<int> N3N_multiplicities;
    // extract N4N multiplicity distribution
    std::vector<int> N4N_multiplicities;
    // extract number of SF reactions, number of IF reactions, number of capture reactions
    int number_of_SF_reactions = 0;
    int number_of_IF_reactions = 0;
    int number_of_N2N_reactions = 0;
    int number_of_N3N_reactions = 0;
    int number_of_N4N_reactions = 0;
    int number_of_B10_captures = 0;
    int number_of_H1_captures = 0;
    int number_of_Cf252_captures = 0;
    // extract the time of B-10 capture
    std::vector<double> B10_capture_times;
    // extract the time of SP source events
    std::vector<double> SF_source_times;
    // extract the time difference between the creation and termination of neutrons
    std::vector<double> neutron_survival_times;
    for (const auto &record : nps_history)
    {
        // extract fission chain data
        FissionChainData fission_chain_data(record);
        for (const auto& neutron_reaction_data : fission_chain_data.reaction_data)
        {
            switch (neutron_reaction_data.reaction_type)
            {
            case SF_REACTION:
                number_of_SF_reactions++;
                SF_multiplicities.push_back(neutron_reaction_data.children_multiplicity);
                break;
            case IF_REACTION:
                number_of_IF_reactions++;
                IF_multiplicities.push_back(neutron_reaction_data.children_multiplicity);
                break;
            case N2N_REACTION:
                number_of_N2N_reactions++;
                N2N_multiplicities.push_back(neutron_reaction_data.children_multiplicity);
                break;
            case N3N_REACTION:
                number_of_N3N_reactions++;
                N3N_multiplicities.push_back(neutron_reaction_data.children_multiplicity);
                break;
            case N4N_REACTION:
                number_of_N4N_reactions++;
                N4N_multiplicities.push_back(neutron_reaction_data.children_multiplicity);
                break;
            case B10_CAPTURE:
                number_of_B10_captures++;
                B10_capture_times.push_back(neutron_reaction_data.reaction_time);
                break;
            case H1_CAPTURE:
                number_of_H1_captures++;
                break;
            case CF252_CAPTURE:
                number_of_Cf252_captures++;
                break;
            default:
                break;
            }
        }
        // append all data in neutron_reaction_data.survival_time to neutron_survival_times
        neutron_survival_times.insert(neutron_survival_times.end(), fission_chain_data.survival_time.begin(), fission_chain_data.survival_time.end());
        // append all data in neutron_reaction_data.SF_event_times to SF_source_times
        SF_source_times.insert(SF_source_times.end(), fission_chain_data.SF_event_times.begin(), fission_chain_data.SF_event_times.end());
    }
    std::vector<int> multiplicity(256, 0);
    std::iota(multiplicity.begin(), multiplicity.end(), 0);
    std::vector<int> SF_multiplicity_counts(256, 0);
    std::vector<int> IF_multiplicity_counts(256, 0);
    std::vector<int> N2N_multiplicity_counts(256, 0);
    std::vector<int> N3N_multiplicity_counts(256, 0);
    std::vector<int> N4N_multiplicity_counts(256, 0);
    for (const auto& multiplicity : SF_multiplicities)
    {
        SF_multiplicity_counts[multiplicity]++;
    }
    for (const auto& multiplicity : IF_multiplicities)
    {
        IF_multiplicity_counts[multiplicity]++;
    }
    for (const auto& multiplicity : N2N_multiplicities)
    {
        N2N_multiplicity_counts[multiplicity]++;
    }
    for (const auto& multiplicity : N3N_multiplicities)
    {
        N3N_multiplicity_counts[multiplicity]++;
    }
    for (const auto& multiplicity : N4N_multiplicities)
    {
        N4N_multiplicity_counts[multiplicity]++;
    }

    // save data to HDF5 file
    std::string output_h5_file_name = "fission_chain_data.h5";
    try {
        // Create a new file using the default property lists.
        HighFive::File file(output_h5_file_name, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Overwrite);

        // Create the dataset
        std::vector<int> number_of_reactions{number_of_SF_reactions, 
                                             number_of_IF_reactions, 
                                             number_of_N2N_reactions, 
                                             number_of_N3N_reactions, 
                                             number_of_N4N_reactions,
                                             number_of_B10_captures, 
                                             number_of_H1_captures, 
                                             number_of_Cf252_captures};
        HighFive::DataSet dataset_number_of_reactions = file.createDataSet<int>("number_of_reactions", HighFive::DataSpace::From(number_of_reactions));
        // write it
        dataset_number_of_reactions.write(number_of_reactions);
        HighFive::DataSet dataset_SF_multiplicity_counts = file.createDataSet<int>("SF_multiplicity_counts", HighFive::DataSpace::From(SF_multiplicity_counts));
        dataset_SF_multiplicity_counts.write(SF_multiplicity_counts);
        HighFive::DataSet dataset_IF_multiplicity_counts = file.createDataSet<int>("IF_multiplicity_counts", HighFive::DataSpace::From(IF_multiplicity_counts));
        dataset_IF_multiplicity_counts.write(IF_multiplicity_counts);
        HighFive::DataSet dataset_N2N_multiplicity_counts = file.createDataSet<int>("N2N_multiplicity_counts", HighFive::DataSpace::From(N2N_multiplicity_counts));
        dataset_N2N_multiplicity_counts.write(N2N_multiplicity_counts);
        HighFive::DataSet dataset_N3N_multiplicity_counts = file.createDataSet<int>("N3N_multiplicity_counts", HighFive::DataSpace::From(N3N_multiplicity_counts));
        dataset_N3N_multiplicity_counts.write(N3N_multiplicity_counts);
        HighFive::DataSet dataset_N4N_multiplicity_counts = file.createDataSet<int>("N4N_multiplicity_counts", HighFive::DataSpace::From(N4N_multiplicity_counts));
        dataset_N4N_multiplicity_counts.write(N4N_multiplicity_counts);
        HighFive::DataSet dataset_B10_capture_times = file.createDataSet<double>("B10_capture_times", HighFive::DataSpace::From(B10_capture_times));
        dataset_B10_capture_times.write(B10_capture_times);
        HighFive::DataSet dataset_neutron_survival_times = file.createDataSet<double>("neutron_survival_times", HighFive::DataSpace::From(neutron_survival_times));
        dataset_neutron_survival_times.write(neutron_survival_times);
        HighFive::DataSet dataset_SF_source_times = file.createDataSet<double>("SF_source_times", HighFive::DataSpace::From(SF_source_times));
        dataset_SF_source_times.write(SF_source_times);
    } catch (HighFive::Exception& err) {
        // catch and print any HDF5 error
        std::cerr << err.what() << std::endl;
    }

    auto endTime = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << "ms" << std::endl; 
    
    return 0;
}