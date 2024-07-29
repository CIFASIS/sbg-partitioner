/*****************************************************************************

 This file is part of SBG Partitioner.

 SBG Partitioner is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 SBG Partitioner is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with SBG Partitioner.  If not, see <http://www.gnu.org/licenses/>.

 ******************************************************************************/

#include <set>
#include <util/logger.hpp>

#include "build_sb_graph.hpp"
#include "dfs_on_sbg.hpp"
#include "partition_graph.hpp"


using namespace std;

using namespace SBG::LIB;
using namespace SBG::Util;

using namespace sbg_partitioner::search;

namespace sbg_partitioner {


PartitionMap make_initial_partition(BaseSBG& graph, unsigned number_of_partitions, PartitionAlgorithm algorithm, bool pre_order)
{
    PartitionMap partitions_set;
    switch (algorithm) {
        case PartitionAlgorithm::DISTRIBUTED:
            initialize_partitioning(graph, number_of_partitions, make_unique<PartitionStrategyDistributive>(number_of_partitions, graph), pre_order);
            break;
        case PartitionAlgorithm::GREEDY:
        default:
            initialize_partitioning(graph, number_of_partitions, make_unique<PartitionStrategyGreedy>(number_of_partitions, graph), pre_order);
    }

    map<unsigned, set<SetPiece>> partitions = partitionate();

    for (auto& [id, set] : partitions) {
        SBG::LIB::UnordSet set_piece;
        for (auto& s : set) {
            if (not s.intervals().empty()) {
                set_piece.emplace(s.intervals().front());
            }
        }
        partitions_set[id] = set_piece;
    }


    SBG_LOG << partitions_set;
    cout << partitions_set << endl;

    return partitions_set;
}


static size_t get_partition_communication(SBG::LIB::BaseSBG& graph, const PartitionMap& partitions)
{
    SBG::LIB::UnordSet s;
    for (auto& [i, _] : partitions) {
        auto ss = get_connectivity_set(graph, partitions, i);
        s = SBG::LIB::cup(ss, s);
        cout << "current connectivity set " << s << ", cardinality " << get_unordset_size(s) << endl;
    }

    size_t size = get_unordset_size(s);

    return size;
}


static void aux_function_best_initial_partition(
    SBG::LIB::BaseSBG& graph,
    PartitionMap& best_initial_partitions,
    size_t& best_communication_set_cardinality,
    unsigned number_of_partitions,
    PartitionAlgorithm partition_algorithm,
    bool pre_order)
{
    PartitionMap temp_intial_partitions =  make_initial_partition(graph, number_of_partitions, PartitionAlgorithm::GREEDY, true);
    size_t temp_partition_comm_size = get_partition_communication(graph, temp_intial_partitions);

    cout  << "PartitionAlgorithm: " << partition_algorithm << ", pre order " << pre_order << " cardinality " << best_communication_set_cardinality << endl;

    if (temp_partition_comm_size < best_communication_set_cardinality) {
        best_initial_partitions = std::move(temp_intial_partitions);
        best_communication_set_cardinality = temp_partition_comm_size;
    }

}


PartitionMap
best_initial_partition(
    SBG::LIB::BaseSBG& graph,
    unsigned number_of_partitions)
{
    PartitionMap best_initial_partitions = make_initial_partition(graph, number_of_partitions, PartitionAlgorithm::DISTRIBUTED, true);
    size_t best_communication_set_cardinality = get_partition_communication(graph, best_initial_partitions);

    cout << "PartitionAlgorithm: " << PartitionAlgorithm::DISTRIBUTED << ", pre order " << true << " cardinality " << best_communication_set_cardinality << endl;

    aux_function_best_initial_partition(graph, best_initial_partitions, best_communication_set_cardinality, number_of_partitions, PartitionAlgorithm::DISTRIBUTED, false);

    aux_function_best_initial_partition(graph, best_initial_partitions, best_communication_set_cardinality, number_of_partitions, PartitionAlgorithm::GREEDY, true);

    aux_function_best_initial_partition(graph, best_initial_partitions, best_communication_set_cardinality, number_of_partitions, PartitionAlgorithm::GREEDY, false);

    cout << "Best is " << best_initial_partitions << " with communication " << best_communication_set_cardinality << endl;

    return best_initial_partitions;
}


static UnordSet get_communication_edges(UnordSet partition, const BasePWMap& map_1, const BasePWMap& map_2)
{
    auto pre_image = preImage(partition, map_1);
    auto image_map_2 = image(pre_image, map_2);
    auto outside_partition = difference(image_map_2, partition);
    auto comm_edges = preImage(outside_partition, map_2);

    return comm_edges;
}


UnordSet get_connectivity_set(
    BaseSBG& graph,
    const PartitionMap& partitions,
    size_t partition_index)
{
    const auto& partition = partitions.at(partition_index);

    UnordSet edges;

    for (const auto& [i, p] : partitions) {
        if (i == partition_index) {
            continue;
        }

        for (size_t i = 0; i < graph.map1().size(); i++) {
            auto map_1 = *(graph.map1().maps().begin() + i);
            auto map_2 = *(graph.map2().maps().begin() + i);
            auto comm_edges_1 = get_communication_edges(partition, map_1, map_2);
            auto comm_edges_2 = get_communication_edges(partition, map_2, map_1);
            auto comm_edges = cup (comm_edges_1, comm_edges_2);
            edges = cup(edges, comm_edges);
        }
    }

    return edges;
}


size_t get_unordset_size(const UnordSet& set)
{
    size_t acc = 0;
    for (auto& set_piece : set.pieces()) {
        for (auto& interval : set_piece.intervals()) {
            acc += (interval.end() - interval.begin() + 1);
        }
    }

    return acc;
}


ostream& operator<<(ostream& os, const PartitionMap& partitions)
{
    for(const auto& [i, p]: partitions) {
        os << i << ", " << p << endl;
    }

    return os;
}

}