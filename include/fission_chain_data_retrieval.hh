#pragma once

#include "fission_chain_tree.hh"


struct ReactionData
{
    int reaction_type=-1;
    double reaction_time=-1.0;
    int children_multiplicity=0;
    // std::array<double, 3> reaction_pos = {0.0, 0.0, 0.0};
    double initial_energy=-1.0;
};


class FissionChainData
{
private:
    void _retrieve_data(const NPSHistory& record)
    {
        // if record.root() has no children, then no neutrons are created in this event, do nothing
        // else:
        if (record.root().size() != 0)
        {
            for (NPSHistoryConstIterator it(record.begin()); it != record.end(); it++)
            {
                ReactionData neutron_reaction_data;
                neutron_reaction_data.reaction_type = it->data().destruction_reaction;
                neutron_reaction_data.reaction_time = it->data().destruction_time;
                neutron_reaction_data.initial_energy = it->data().creation_energy;
                neutron_reaction_data.children_multiplicity = it->size();
                if(!it->is_root())
                {
                    survival_time.push_back(it->data().destruction_time - it->data().creation_time);
                }
                else
                {
                    SF_event_times.push_back(it->data().destruction_time);
                }
                reaction_data.push_back(neutron_reaction_data);
            }
        }
    }
public:
    std::vector<ReactionData> reaction_data;
    std::vector<double> survival_time;
    std::vector<double> SF_event_times;
    FissionChainData(const NPSHistory& record)
    {
        _retrieve_data(record);
    }
};