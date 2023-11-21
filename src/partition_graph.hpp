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

class PartitionGraph {

public: 
    PartitionGraph(SBG::LIB::CanonSBG& _graph, int number_of_partitions);

    const SBG::LIB::CanonSBG& graph() const;

    void set_partition(size_t interval_index, size_t partition_number);

    std::map<SBG::LIB::SetPiece, int> partitions() const;

    std::unordered_set<size_t> get_connectivity_set(size_t edge_index) const;

private:
    SBG::LIB::CanonSBG _graph;
    std::map<SBG::LIB::SetPiece, int> _partitions;

    void make_initial_partition(int number_of_partitions);
};

size_t connectivity_set_cardinality(const PartitionGraph& pgraph, size_t edge_index);

std::ostream& operator<<(std::ostream& os, const PartitionGraph pgraph);

}