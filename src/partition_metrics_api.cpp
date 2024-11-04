/**
 This file is part of Set--Based Graph Library.

 SBG Library is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SBG Library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with SBG Library.  If not, see <http://www.gnu.org/licenses/>.

 ******************************************************************************/

#include <algorithm>
#include <fstream>
#include <set>

#include <sbg/sbg.hpp>

#include "build_sb_graph.hpp"
#include "partition_metrics_api.hpp"
#include "weighted_sb_graph.hpp"

using namespace std;

using namespace SBG::LIB;
using namespace SBG::Util;

namespace sbg_partitioner{

namespace metrics {

namespace {

UnordSet get_edge_cut(
    const UnordSet& partition_a,
    const UnordSet& partition_b,
    const MapSet<SBG::LIB::UnordSet>& maps_1,
    const MapSet<SBG::LIB::UnordSet>& maps_2)
{
    UnordSet external_communication;
    for (size_t i = 0; i < maps_1.size(); i++) {
        auto map_1 = *(maps_1.begin() + i);
        auto map_2 = *(maps_2.begin() + i);
        auto d = preImage(partition_a, map_1);
        auto im = image(d, map_2);
        auto ec_nodes = intersection(partition_b, im);
        auto _external_communication = preImage(ec_nodes, map_2);
        external_communication = cup(_external_communication, external_communication);
    }

    return external_communication;
}


int communication_volume_one_dim(const PartitionMap& partitions, const SetPiece& v, const WeightedSBGraph& sb_graph, unsigned i)
{
    int acc = 0;
    for (unsigned value = v.intervals()[0].begin(); value <= v.intervals()[0].end(); value += v.intervals()[0].step()) {
        int local_total_vol = 0;
        Interval set_piece = Interval(value, v.intervals()[0].step(), value);
        const UnordSet adjacents = get_adjacents(sb_graph, set_piece);

        for (unsigned j = 0; j < partitions.size(); j++) {
            if (i == j) {
                continue;
            }

            const auto& p_set = partitions.at(j);

            if (not isEmpty(intersection(adjacents, p_set))) {
                local_total_vol += 1;
            }
        }

        acc += local_total_vol;
    }

    return acc;
}


int communication_volume_two_dim(const PartitionMap& partitions, const SetPiece& v, const WeightedSBGraph& sb_graph, unsigned i)
{
    int acc = 0;
    for (unsigned v_0 = v.intervals()[0].begin(); v_0 <= v.intervals()[0].end(); v_0 += v.intervals()[0].step()) {
        int local_total_vol = 0;
        Interval interval_0 = Interval(v_0, v.intervals()[0].step(), v_0);
        for (unsigned v_1 = v.intervals()[1].begin(); v_1 <= v.intervals()[1].end(); v_1 += v.intervals()[1].step()) {
            Interval interval_1 = Interval(v_1, v.intervals()[1].step(), v_1);

            SetPiece set_piece;
            set_piece.emplaceBack(interval_0);
            set_piece.emplaceBack(interval_1);

            const UnordSet adjacents = get_adjacents(sb_graph, set_piece);

            for (unsigned j = 0; j < partitions.size(); j++) {
                if (i == j) {
                    continue;
                }

                const auto& p_set = partitions.at(j);
                if (not isEmpty(intersection(adjacents, p_set))) {
                    local_total_vol += 1;
                }
            }

            acc += local_total_vol;
        }
    }

    return acc;
}


void write_node_by_partition(const PartitionMap& partitions, const WeightedSBGraph& sb_graph)
{
    vector<SetPiece> nodes;
    nodes.reserve(sb_graph.V().size());
    for (auto v : sb_graph.V()) {
        nodes.push_back(v);
    }

    // SORT BIDIMENSIONAL:
    if (sb_graph.V()[0].size() == 2) {
        auto f_sort = [] (const auto& a, const auto& b) {
            if (a.intervals()[0].end() < b.intervals()[0].end()) {
                return true;
            } else if (b.intervals()[0].end() < a.intervals()[0].end()) {
                return false;
            }

            return a.intervals()[1].end() < b.intervals()[1].end();
        };

        sort(nodes.begin(), nodes.end(), f_sort);

        // expand it and write it
        vector<unsigned> partition_by_node;
        ofstream output_file("output.txt");
        for (const auto& n : nodes) {
            for (unsigned v_0 = n.intervals()[0].begin(); v_0 <= n.intervals()[0].end(); v_0++) {
                for (unsigned v_1 = n.intervals()[1].begin(); v_1 <= n.intervals()[1].end(); v_1++) {
                    SetPiece set_piece;
                    set_piece.emplaceBack(Interval(v_0, n.intervals()[0].step(), v_0));
                    set_piece.emplaceBack(Interval(v_1, n.intervals()[1].step(), v_1));
                    for (const auto& [i, p] : partitions) {
                        if (not isEmpty(intersection(set_piece, p))) {
                            partition_by_node.push_back(i);
                            output_file << to_string(i) << endl;
                            break;
                        }
                    }
                }
            }
        }
    }

    // SORT ONE-DIMENSIONAL:
    if (sb_graph.V()[0].size() == 1) {
        auto f_sort = [] (const auto& a, const auto& b) {
            return a.intervals()[0].end() < b.intervals()[0].end();
        };

        sort(nodes.begin(), nodes.end(), f_sort);
        cout << "sorted nodes? ";
        for (const auto& v : nodes) {
            cout << v << " ";
        }
        cout << endl;
    }
}

}

int edge_cut(const PartitionMap& partitions, const WeightedSBGraph& sb_graph)
{
    UnordSet ec;
    const auto& maps_1 = sb_graph.map1().maps();
    const auto& maps_2 = sb_graph.map2().maps();;
    for (size_t i = 0; i < partitions.size(); i++) {
        const UnordSet& partition_1 = partitions.at(i);
        for (size_t j = i + 1; j < partitions.size(); j++) {
            const UnordSet& partition_2 = partitions.at(j);
            ec = cup(get_edge_cut(partition_1, partition_2, maps_1, maps_2), ec);
            ec = cup(get_edge_cut(partition_1, partition_2, maps_2, maps_1), ec);
        }
    }

    int weight = get_edge_set_cost(ec, sb_graph.get_edge_costs());

    return weight;
}


pair<int, int> communication_volume(const PartitionMap& partitions, const WeightedSBGraph& sb_graph)
{
    int comm_vol = 0;
    int max_comm_vol = 0;
    for (unsigned i = 0; i < partitions.size(); i++) {
        int communication_volume_partition = 0;
        const auto& p = partitions.at(i);

        for (const auto& v : p.pieces()) {
            if (v.size() == 1) {
                communication_volume_partition += communication_volume_one_dim(partitions, v, sb_graph, i);
            } else if (v.size() == 2) {
                communication_volume_partition += communication_volume_two_dim(partitions, v, sb_graph, i);
            } else {
                throw 2;
            }
        }

        comm_vol += communication_volume_partition;
        max_comm_vol = max(communication_volume_partition, max_comm_vol);
    }

    return {comm_vol, max_comm_vol};
}


float maximum_imbalance(const PartitionMap& partitions, const WeightedSBGraph& sb_graph)
{
    unsigned number_of_nodes = get_node_size(sb_graph.V(), sb_graph.get_node_weights());
    float expected_imb = number_of_nodes / partitions.size();

    float max_imbalance = 0.;
    for (const auto& [i, p] : partitions) {
        unsigned size_of_p = get_node_size(p, sb_graph.get_node_weights());
        float imbalance_p = abs(expected_imb - float(size_of_p)) / expected_imb;
        max_imbalance = max(max_imbalance, imbalance_p);
    }

    write_node_by_partition(partitions, sb_graph);

    return max_imbalance;
}

}

}
