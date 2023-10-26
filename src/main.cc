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

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

#include "parser.hh"
int main(int argc, char** argv)
{
    auto startTime = std::chrono::high_resolution_clock::now();

    // const std::string ptracFilePath(argv[1]);
    const std::string ptracFilePath("/media/ming/Elements/multiplicity_test/BCS_NMC_1500/cf252_sphere/ptrac/v11/ptrac");
    MCNPPTRACBinary ptracFile(ptracFilePath);
    const long maxNum(1e9);
    // std::vector<ParticleHistory> spontaneous_fission_neutrons;
    // std::vector<ParticleHistory> induced_fission_neutrons;
    std::vector<long> SF_neutron_nps;
    std::vector<int>  SF_neutron_creation_eventID;
    std::vector<double> SF_neutron_creation_pos_x;
    std::vector<double> SF_neutron_creation_pos_y;
    std::vector<double> SF_neutron_creation_pos_z;
    std::vector<double> SF_neutron_creation_energy;
    std::vector<double> SF_neutron_creation_time;
    std::vector<int>  SF_neutron_termination_eventID;
    std::vector<double> SF_neutron_termination_pos_x;
    std::vector<double> SF_neutron_termination_pos_y;
    std::vector<double> SF_neutron_termination_pos_z;
    std::vector<double> SF_neutron_termination_energy;
    std::vector<double> SF_neutron_termination_time;

    std::vector<long> IF_neutron_nps;
    std::vector<int>  IF_neutron_creation_eventID;
    std::vector<double> IF_neutron_creation_pos_x;
    std::vector<double> IF_neutron_creation_pos_y;
    std::vector<double> IF_neutron_creation_pos_z;
    std::vector<double> IF_neutron_creation_energy;
    std::vector<double> IF_neutron_creation_time;
    std::vector<int>  IF_neutron_termination_eventID;
    std::vector<double> IF_neutron_termination_pos_x;
    std::vector<double> IF_neutron_termination_pos_y;
    std::vector<double> IF_neutron_termination_pos_z;
    std::vector<double> IF_neutron_termination_energy;
    std::vector<double> IF_neutron_termination_time;

    bool flag = true;
    while (ptracFile.readNextNPS(maxNum))
    {
        
        if (ptracFile.getNPSRead() % 10000 == 0)
        {
            std::cout << "NPS = " << ptracFile.getNPSRead() << '\n';
        }
        const NPSHistory record = ptracFile.getNPSHistory();
        for (auto && neutron_history : record)
        {
            const auto& creation_event = neutron_history[0];
            const auto& termination_event = neutron_history[1];
            // if eventID is 1000 or 2000, then it is a SP fission neutron
            if (creation_event.eventID == 1000 || creation_event.eventID == 2000)
            {
                // spontaneous_fission_neutrons.push_back(neutron_history);
                SF_neutron_nps.push_back(creation_event.nps);
                SF_neutron_creation_eventID.push_back(creation_event.eventID);
                SF_neutron_creation_pos_x.push_back(creation_event.pos[0]);
                SF_neutron_creation_pos_y.push_back(creation_event.pos[1]);
                SF_neutron_creation_pos_z.push_back(creation_event.pos[2]);
                SF_neutron_creation_energy.push_back(creation_event.energy);
                SF_neutron_creation_time.push_back(creation_event.time);

                SF_neutron_termination_eventID.push_back(termination_event.eventID);
                SF_neutron_termination_pos_x.push_back(termination_event.pos[0]);
                SF_neutron_termination_pos_y.push_back(termination_event.pos[1]);
                SF_neutron_termination_pos_z.push_back(termination_event.pos[2]);
                SF_neutron_termination_energy.push_back(termination_event.energy);
                SF_neutron_termination_time.push_back(termination_event.time);  
            }
            // if eventID is 2007, then it is a induced fission neutron
            else if (creation_event.eventID == 2007)
            {
                // induced_fission_neutrons.push_back(neutron_history);
                IF_neutron_nps.push_back(creation_event.nps);
                IF_neutron_creation_eventID.push_back(creation_event.eventID);
                IF_neutron_creation_pos_x.push_back(creation_event.pos[0]);
                IF_neutron_creation_pos_y.push_back(creation_event.pos[1]);
                IF_neutron_creation_pos_z.push_back(creation_event.pos[2]);
                IF_neutron_creation_energy.push_back(creation_event.energy);
                IF_neutron_creation_time.push_back(creation_event.time);

                IF_neutron_termination_eventID.push_back(termination_event.eventID);
                IF_neutron_termination_pos_x.push_back(termination_event.pos[0]);
                IF_neutron_termination_pos_y.push_back(termination_event.pos[1]);
                IF_neutron_termination_pos_z.push_back(termination_event.pos[2]);
                IF_neutron_termination_energy.push_back(termination_event.energy);
                IF_neutron_termination_time.push_back(termination_event.time);
            }
        }
    }

    // write fission neutrons to HDF5 file
    std::string output_h5_file_name = "fission_neutron_histories.h5";
    try {
        // Create a new file using the default property lists.
        HighFive::File file(output_h5_file_name, HighFive::File::ReadWrite | HighFive::File::Create | HighFive::File::Overwrite);

        // Create the dataset
        HighFive::DataSet SF_nps_dataset = file.createDataSet<long>("SF/nps", HighFive::DataSpace::From(SF_neutron_nps));
        SF_nps_dataset.write(SF_neutron_nps);
        HighFive::DataSet SF_creation_eventID_dataset = file.createDataSet<int>("SF/creation_eventID", HighFive::DataSpace::From(SF_neutron_creation_eventID));
        SF_creation_eventID_dataset.write(SF_neutron_creation_eventID);
        HighFive::DataSet SF_creation_pos_x_dataset = file.createDataSet<double>("SF/creation_pos_x", HighFive::DataSpace::From(SF_neutron_creation_pos_x));
        SF_creation_pos_x_dataset.write(SF_neutron_creation_pos_x);
        HighFive::DataSet SF_creation_pos_y_dataset = file.createDataSet<double>("SF/creation_pos_y", HighFive::DataSpace::From(SF_neutron_creation_pos_y));
        SF_creation_pos_y_dataset.write(SF_neutron_creation_pos_y);
        HighFive::DataSet SF_creation_pos_z_dataset = file.createDataSet<double>("SF/creation_pos_z", HighFive::DataSpace::From(SF_neutron_creation_pos_z));
        SF_creation_pos_z_dataset.write(SF_neutron_creation_pos_z);
        HighFive::DataSet SF_creation_energy_dataset = file.createDataSet<double>("SF/creation_energy", HighFive::DataSpace::From(SF_neutron_creation_energy));
        SF_creation_energy_dataset.write(SF_neutron_creation_energy);
        HighFive::DataSet SF_creation_time_dataset = file.createDataSet<double>("SF/creation_time", HighFive::DataSpace::From(SF_neutron_creation_time));
        SF_creation_time_dataset.write(SF_neutron_creation_time);

        HighFive::DataSet SF_termination_eventID_dataset = file.createDataSet<int>("SF/termination_eventID", HighFive::DataSpace::From(SF_neutron_termination_eventID));
        SF_termination_eventID_dataset.write(SF_neutron_termination_eventID);
        HighFive::DataSet SF_termination_pos_x_dataset = file.createDataSet<double>("SF/termination_pos_x", HighFive::DataSpace::From(SF_neutron_termination_pos_x));
        SF_termination_pos_x_dataset.write(SF_neutron_termination_pos_x);
        HighFive::DataSet SF_termination_pos_y_dataset = file.createDataSet<double>("SF/termination_pos_y", HighFive::DataSpace::From(SF_neutron_termination_pos_y));
        SF_termination_pos_y_dataset.write(SF_neutron_termination_pos_y);
        HighFive::DataSet SF_termination_pos_z_dataset = file.createDataSet<double>("SF/termination_pos_z", HighFive::DataSpace::From(SF_neutron_termination_pos_z));
        SF_termination_pos_z_dataset.write(SF_neutron_termination_pos_z);
        HighFive::DataSet SF_termination_energy_dataset = file.createDataSet<double>("SF/termination_energy", HighFive::DataSpace::From(SF_neutron_termination_energy));
        SF_termination_energy_dataset.write(SF_neutron_termination_energy);
        HighFive::DataSet SF_termination_time_dataset = file.createDataSet<double>("SF/termination_time", HighFive::DataSpace::From(SF_neutron_termination_time));
        SF_termination_time_dataset.write(SF_neutron_termination_time);

        HighFive::DataSet IF_nps_dataset = file.createDataSet<long>("IF/nps", HighFive::DataSpace::From(IF_neutron_nps));
        IF_nps_dataset.write(IF_neutron_nps);
        HighFive::DataSet IF_creation_eventID_dataset = file.createDataSet<int>("IF/creation_eventID", HighFive::DataSpace::From(IF_neutron_creation_eventID));
        IF_creation_eventID_dataset.write(IF_neutron_creation_eventID);
        HighFive::DataSet IF_creation_pos_x_dataset = file.createDataSet<double>("IF/creation_pos_x", HighFive::DataSpace::From(IF_neutron_creation_pos_x));
        IF_creation_pos_x_dataset.write(IF_neutron_creation_pos_x);
        HighFive::DataSet IF_creation_pos_y_dataset = file.createDataSet<double>("IF/creation_pos_y", HighFive::DataSpace::From(IF_neutron_creation_pos_y));
        IF_creation_pos_y_dataset.write(IF_neutron_creation_pos_y);
        HighFive::DataSet IF_creation_pos_z_dataset = file.createDataSet<double>("IF/creation_pos_z", HighFive::DataSpace::From(IF_neutron_creation_pos_z));
        IF_creation_pos_z_dataset.write(IF_neutron_creation_pos_z);
        HighFive::DataSet IF_creation_energy_dataset = file.createDataSet<double>("IF/creation_energy", HighFive::DataSpace::From(IF_neutron_creation_energy));
        IF_creation_energy_dataset.write(IF_neutron_creation_energy);
        HighFive::DataSet IF_creation_time_dataset = file.createDataSet<double>("IF/creation_time", HighFive::DataSpace::From(IF_neutron_creation_time));
        IF_creation_time_dataset.write(IF_neutron_creation_time);

        HighFive::DataSet IF_termination_eventID_dataset = file.createDataSet<int>("IF/termination_eventID", HighFive::DataSpace::From(IF_neutron_termination_eventID));
        IF_termination_eventID_dataset.write(IF_neutron_termination_eventID);
        HighFive::DataSet IF_termination_pos_x_dataset = file.createDataSet<double>("IF/termination_pos_x", HighFive::DataSpace::From(IF_neutron_termination_pos_x));
        IF_termination_pos_x_dataset.write(IF_neutron_termination_pos_x);
        HighFive::DataSet IF_termination_pos_y_dataset = file.createDataSet<double>("IF/termination_pos_y", HighFive::DataSpace::From(IF_neutron_termination_pos_y));
        IF_termination_pos_y_dataset.write(IF_neutron_termination_pos_y);
        HighFive::DataSet IF_termination_pos_z_dataset = file.createDataSet<double>("IF/termination_pos_z", HighFive::DataSpace::From(IF_neutron_termination_pos_z));
        IF_termination_pos_z_dataset.write(IF_neutron_termination_pos_z);
        HighFive::DataSet IF_termination_energy_dataset = file.createDataSet<double>("IF/termination_energy", HighFive::DataSpace::From(IF_neutron_termination_energy));
        IF_termination_energy_dataset.write(IF_neutron_termination_energy);
        HighFive::DataSet IF_termination_time_dataset = file.createDataSet<double>("IF/termination_time", HighFive::DataSpace::From(IF_neutron_termination_time));
        IF_termination_time_dataset.write(IF_neutron_termination_time);

    } catch (HighFive::Exception& err) {
        // catch and print any HDF5 error
        std::cerr << err.what() << std::endl;
    }
    auto endTime = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count() << "ms" << std::endl; 
    
    return 0;
}