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

#pragma once

#include <set>


#include "partition_graph.hpp"
#include "weighted_sb_graph.hpp"


namespace sbg_partitioner {

typedef std::map<unsigned, int> ImbalanceMap;


struct GainObjectImbalance {
    size_t i;
    size_t j;
    int gain;
    size_t size_i;
    size_t size_j;
};

std::ostream& operator<<(std::ostream& os, const GainObjectImbalance& gain);

using GainObjectImbalanceComparator = GainObjectComparatorTemplate<GainObjectImbalance>;

using CostMatrixImbalance = std::set<GainObjectImbalance, GainObjectImbalanceComparator>;

std::ostream& operator<<(std::ostream& os, const CostMatrixImbalance& cost_matrix);

std::pair<unsigned, unsigned>
compute_lmin_lmax(const WeightedSBGraph& graph, unsigned number_of_partitions, const float imbalance_epsilon);


void compute_partition_imbalance(
    unsigned i, unsigned j,
    SBG::LIB::UnordSet& partition_a,
    SBG::LIB::UnordSet& partition_b,
    const WeightedSBGraph& graph,
    const NodeWeight& node_weight,
    unsigned LMin,
    unsigned LMax,
    CostMatrixImbalance& cost_matrix
);


CostMatrixImbalance generate_gain_matrix(
    const WeightedSBGraph& graph,
    const NodeWeight& node_weight,
    SBG::LIB::UnordSet& partition_a,
    SBG::LIB::UnordSet& partition_b,
    unsigned LMin,
    unsigned LMax);


int kl_sbg_imbalance(
    const WeightedSBGraph& graph,
    SBG::LIB::UnordSet& partition_a,
    SBG::LIB::UnordSet& partition_b,
    unsigned LMin,
    unsigned LMax);


KLBipartResult kl_sbg_bipart_imbalance(
    const WeightedSBGraph& graph,
    SBG::LIB::UnordSet& partition_a,
    SBG::LIB::UnordSet& partition_b,
    unsigned LMin,
    unsigned LMax);


void kl_sbg_imbalance_partitioner(const WeightedSBGraph& graph, PartitionMap& partitions, const float imbalance_epsilon);

}