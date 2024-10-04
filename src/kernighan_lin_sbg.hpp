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

#include <set>

#include "build_sb_graph.hpp"
#include "partition_graph.hpp"

namespace sbg_partitioner
{



struct GainObject {
    size_t i;
    size_t j;
    int gain;
    size_t size;
};





using GainObjectComparator = GainObjectComparatorTemplate<GainObject>;

using CostMatrix = std::set<GainObject, GainObjectComparator>;





void kl_sbg_partitioner(const WeightedSBGraph& graph, PartitionMap& partitions);

KLBipartResult kl_sbg_bipart(const WeightedSBGraph& graph, SBG::LIB::UnordSet& partition_a, SBG::LIB::UnordSet& partition_b);

int kl_sbg(const WeightedSBGraph& graph, SBG::LIB::UnordSet& partition_a, SBG::LIB::UnordSet& partition_b);

void update_sum(
    int& par_sum,
    int g,
    int& max_par_sum,
    std::pair<SBG::LIB::UnordSet, SBG::LIB::UnordSet>& max_par_sum_set,
    SBG::LIB::UnordSet& a_v,
    SBG::LIB::UnordSet& b_v);



std::ostream& operator<<(std::ostream& os, const GainObject& gain);

std::ostream& operator<<(std::ostream& os, const CostMatrix& cost_matrix);




}; // namespace sbg_partitioner
