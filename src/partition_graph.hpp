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

#include <iostream>
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


/// This class represents a partioned set based graph.
/// The main purpose of PartitionGraph is to divide the
/// calculation of different nodes of a dynamic model
/// in order to simulate it in parallel.
class PartitionGraph {

public:
    PartitionGraph(SBG::LIB::CanonSBG& _graph, unsigned number_of_partitions,
                   PartitionAlgorithm algorithm, bool pre_order);

    /// Returns a reference of the set based graph to be partioned.
    const SBG::LIB::CanonSBG& graph() const;

    /// Returns the current partion as a map where each interval is associated
    /// to a partition number.
    std::map<unsigned, SBG::LIB::OrdSet> partitions() const;

    /// Returns the connectivity set of a set of edges contained in map1 and map2 of
    /// the graph (I mean, edges in CanonSBG::map1()[edge_index] and CanonSBG::map2()[edge_index]).
    /// So that, we consider the graph as an undirected graph.
    std::unordered_set<size_t> get_connectivity_set(size_t edge_index) const;

    /// Takes a node and its partition number and returns the internal cost according
    /// to the definition in "An Efficient Heuristic Procedure for Partitioning Graphs".
    /// Note: if node does not fully belong to "partition", it will trigger an exception.
    unsigned int internal_cost(const SBG::LIB::SetPiece& node, unsigned partition);

    /// Takes a node and its partition number and returns the external cost according
    /// to the definition in "An Efficient Heuristic Procedure for Partitioning Graphs"
    /// Note: if node totally or partially belongs to "partition", it will trigger an exception.
    unsigned int external_cost(const SBG::LIB::SetPiece& node, unsigned partition);

    /// Gets you the gain (the reduction in cost) of exchanging nodes A and B.
    /// @param node_a pair with partition number and a subset of that partition
    /// @param node_b analog node_a
    /// @return the gain of exchanging them
    int gain(const std::pair<unsigned, SBG::LIB::SetPiece> node_a, const std::pair<unsigned, SBG::LIB::SetPiece> node_b);

private:
    SBG::LIB::CanonSBG _graph;
    std::map<unsigned, SBG::LIB::OrdSet> _partitions;

    /// It creates an initial partition.
    void make_initial_partition(unsigned number_of_partitions, PartitionAlgorithm algorithm, bool pre_order);
};

/// This function returns the cardinality of a connectivity set.
/// @param pgraph graph partition
/// @param edge_index index of the edges to get the connectivity set.
/// @return the cardinality of the connectivity set
size_t connectivity_set_cardinality(const PartitionGraph& pgraph, size_t edge_index);


std::ostream& operator<<(std::ostream& os, const PartitionGraph& pgraph);

}