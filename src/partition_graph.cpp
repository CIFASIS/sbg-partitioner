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

#include <rapidjson/document.h>
#include <rapidjson/filewritestream.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/writer.h>
#include <set>
#include <util/logger.hpp>

#include "build_sb_graph.hpp"
#include "dfs_on_sbg.hpp"
#include "partition_graph.hpp"


#define PARTITION_SANITY_CHECK 1


using namespace std;

using namespace SBG::LIB;
using namespace SBG::Util;

using namespace sbg_partitioner::search;

namespace sbg_partitioner {


PartitionMap make_initial_partition(WeightedSBGraph& graph, unsigned number_of_partitions, PartitionAlgorithm algorithm, bool pre_order)
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
            SBG::LIB::SetPiece intervals;
            if (not s.intervals().empty()) {
                for (size_t i = 0; i < s.intervals().size(); i++) {
                    Interval interv = s.intervals()[i];
                    intervals.emplaceBack(interv);
                }
            }
            set_piece.emplace(intervals);
        }
        partitions_set[id] = set_piece;
    }

    sanity_check(graph, partitions_set, number_of_partitions);

    SBG_LOG << partitions_set;
    cout << partitions_set << endl;

    return partitions_set;
}


static size_t get_partition_communication(WeightedSBGraph& graph, const PartitionMap& partitions)
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
    WeightedSBGraph& graph,
    PartitionMap& best_initial_partitions,
    size_t& best_communication_set_cardinality,
    unsigned number_of_partitions,
    PartitionAlgorithm partition_algorithm,
    bool pre_order)
{
    PartitionMap temp_intial_partitions =  make_initial_partition(graph, number_of_partitions, partition_algorithm, pre_order);
    size_t temp_partition_comm_size = get_partition_communication(graph, temp_intial_partitions);

    cout  << "PartitionAlgorithm: " << partition_algorithm << ", pre order " << pre_order << " cardinality " << best_communication_set_cardinality << endl;

    if (temp_partition_comm_size < best_communication_set_cardinality) {
        best_initial_partitions = std::move(temp_intial_partitions);
        best_communication_set_cardinality = temp_partition_comm_size;
    }
}


PartitionMap
best_initial_partition(
    WeightedSBGraph& graph,
    unsigned number_of_partitions)
{
    constexpr bool pre_order = true;
    PartitionMap best_initial_partitions = make_initial_partition(graph, number_of_partitions, PartitionAlgorithm::DISTRIBUTED, pre_order);
    size_t best_communication_set_cardinality = get_partition_communication(graph, best_initial_partitions);

    cout << "PartitionAlgorithm: " << PartitionAlgorithm::DISTRIBUTED << ", pre order " << pre_order << " cardinality " << best_communication_set_cardinality << endl;

    aux_function_best_initial_partition(graph, best_initial_partitions, best_communication_set_cardinality, number_of_partitions, PartitionAlgorithm::DISTRIBUTED, not pre_order);

    aux_function_best_initial_partition(graph, best_initial_partitions, best_communication_set_cardinality, number_of_partitions, PartitionAlgorithm::GREEDY, pre_order);

    aux_function_best_initial_partition(graph, best_initial_partitions, best_communication_set_cardinality, number_of_partitions, PartitionAlgorithm::GREEDY, not pre_order);

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


void sanity_check(const WeightedSBGraph &graph, PartitionMap& partitions_set, unsigned number_of_partitions)
{
# if PARTITION_SANITY_CHECK
    // This is just a sanity check
    UnordSet nodes_to_check;
    for (unsigned i = 0; i < number_of_partitions; i++) {
        nodes_to_check = cup(nodes_to_check, partitions_set[i]);
    }
    UnordSet diff_1 = difference(graph.V(), nodes_to_check);
    UnordSet diff_2 = difference(nodes_to_check, graph.V());
    assert(get_node_size(diff_1) == 0 and "The intial partition has less elements than the graph");
    assert(get_node_size(diff_2) == 0 and "The intial partition has more elements than the graph");
    for (unsigned i = 0; i < number_of_partitions; i++) {
        for (unsigned j = i + 1; j < number_of_partitions; j++) {
            auto p_1 = partitions_set[i];
            auto p_2 = partitions_set[j];
            stringstream error_msg;
            error_msg << "Intersection between " << i << " and " << j << " is not empty." << endl;
            assert(intersection(p_1, p_2).pieces().empty() and error_msg.str().c_str());
        }
    }
#endif //PARTITION_SANITY_CHECK
}


void write_output(const string filename, const PartitionMap& partition_map)
{
    rapidjson::Document json_doc;
    rapidjson::Document::AllocatorType& allocator = json_doc.GetAllocator();
    json_doc.SetObject();

    rapidjson::Value obj_partitions(rapidjson::kArrayType);
    for (size_t i = 0; i < partition_map.size(); i++) {
        rapidjson::Value obj_partition(rapidjson::kArrayType);
        const auto& partition = partition_map.at(i);
        for (const SetPiece& set_piece : partition.pieces()) {
            rapidjson::Value obj_intervals(rapidjson::kArrayType);
            obj_intervals.SetArray();
            for (const Interval& interval : set_piece.intervals()) {
                rapidjson::Value obj_interval(rapidjson::kArrayType);
                rapidjson::Value begin(rapidjson::kNumberType);
                begin.SetUint(interval.begin());
                obj_interval.PushBack(begin, allocator);
                rapidjson::Value end(rapidjson::kNumberType);
                end.SetUint(interval.end());
                obj_interval.PushBack(end, allocator);

                obj_intervals.PushBack(obj_interval, allocator);
            }

            obj_partition.PushBack(obj_intervals, allocator);
        }

        rapidjson::Value obj_nodes(rapidjson::kObjectType);
        obj_nodes.AddMember("nodes", obj_partition, allocator);

        obj_partitions.PushBack(obj_nodes, allocator);
    }

    json_doc.AddMember("partitions", obj_partitions, allocator);

    // Write the JSON data to the file
    rapidjson::StringBuffer s;
    rapidjson::Writer<rapidjson::StringBuffer> writer(s);
    json_doc.Accept(writer);
    std::cout << s.GetString() << std::endl;

    std::ofstream of (filename);
    of << s.GetString();
}


ostream& operator<<(ostream& os, const PartitionMap& partitions)
{
    for(const auto& [i, p]: partitions) {
        os << i << ", " << p << endl;
    }

    return os;
}

}