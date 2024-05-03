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


namespace sbg_partitioner {

enum PartitionAlgorithm
{
    GREEDY = 0,
    DISTRIBUTED = 1
};

// I wish this was a separate function, not part of PartitionGraph but there were a lot of
// compile problems if partitions map object is created locally and OrdSet objects are added.
std::map<unsigned, SBG::LIB::OrdSet>
make_initial_partition(
    SBG::LIB::CanonSBG& graph,
    unsigned number_of_partitions,
    PartitionAlgorithm algorithm,
    bool pre_order);


/// Returns the connectivity set of a set of edges contained in map1 and map2 of
/// the graph (I mean, edges in CanonSBG::map1()[edge_index] and CanonSBG::map2()[edge_index]).
/// So that, we consider the graph as an undirected graph.
std::unordered_set<size_t> get_connectivity_set(
    SBG::LIB::CanonSBG& graph,
    std::map<unsigned, SBG::LIB::OrdSet>& partitions,
    size_t edge_index);


/// This function returns the cardinality of a connectivity set.
/// @param pgraph graph partition
/// @param edge_index index of the edges to get the connectivity set.
/// @return the cardinality of the connectivity set
size_t connectivity_set_cardinality(
    SBG::LIB::CanonSBG& graph,
    std::map<unsigned, SBG::LIB::OrdSet> partitions,
    size_t edge_index);

}