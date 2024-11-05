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

#pragma once

#include <map>
#include <unordered_set>

#include <sbg/interval.hpp>
#include <sbg/sbg.hpp>

#include "weighted_sb_graph.hpp"


namespace sbg_partitioner {

typedef std::map<unsigned, SBG::LIB::OrdSet> PartitionMap;

enum PartitionAlgorithm
{
    GREEDY = 0,
    DISTRIBUTED = 1
};

// I wish this was a separate function, not part of PartitionGraph but there were a lot of
// compile problems if partitions map object is created locally and OrdSet objects are added.
std::vector<PartitionMap>
make_initial_partitions(
    WeightedSBGraph& graph,
    unsigned number_of_partitions);


PartitionMap
best_initial_partition(
    WeightedSBGraph& graph,
    unsigned number_of_partitions);


/// Returns the connectivity set of a set of edges contained in map1 and map2 of
/// the graph (I mean, edges in CanonSBG::map1()[edge_index] and CanonSBG::map2()[edge_index]).
/// So that, we consider the graph as an undirected graph.
SBG::LIB::OrdSet get_connectivity_set(
    SBG::LIB::CanonSBG& graph,
    const PartitionMap& partitions,
    size_t edge_index);


/// This function returns the cardinality of a OrdSet.
size_t get_OrdSet_size(const SBG::LIB::OrdSet& set);


std::string get_output(const PartitionMap& partition_map);


void sanity_check(const WeightedSBGraph& graph, PartitionMap& partitions_set, unsigned number_of_partitions);


std::ostream& operator<<(std::ostream& os, const PartitionMap& partitions);

}