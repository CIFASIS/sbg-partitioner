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

constexpr bool multithreading_enabled = false;

struct GainObject {
    size_t i;
    size_t j;
    int gain;
    size_t size;
};

struct KLBipartResult {
    SBG::LIB::UnordSet A;
    SBG::LIB::UnordSet B;
    int gain;
};


using ec_ic = std::pair<SBG::LIB::UnordPWMDInter , SBG::LIB::UnordPWMDInter>;

using GainObjectComparator = GainObjectComparatorTemplate<GainObject>;

using CostMatrix = std::set<GainObject, GainObjectComparator>;


struct kl_sbg_partitioner_result
{
    size_t i;
    size_t j;
    int gain;
    SBG::LIB::UnordSet A;
    SBG::LIB::UnordSet B;
};


void kl_sbg_partitioner(const WeightedSBGraph& graph, PartitionMap& partitions);

KLBipartResult kl_sbg_bipart(const WeightedSBGraph& graph, SBG::LIB::UnordSet& partition_a, SBG::LIB::UnordSet& partition_b);

int kl_sbg(const WeightedSBGraph& graph, SBG::LIB::UnordSet& partition_a, SBG::LIB::UnordSet& partition_b);

ec_ic compute_EC_IC(
    const SBG::LIB::UnordSet& partition,
    const SBG::LIB::UnordSet& nodes,
    const SBG::LIB::UnordSet& partition_2,
    const SBG::LIB::BasePWMap& departure_map,
    const SBG::LIB::BasePWMap& arrival_map);


void update_sum(
    int& par_sum,
    int g,
    int& max_par_sum,
    std::pair<SBG::LIB::UnordSet, SBG::LIB::UnordSet>& max_par_sum_set,
    SBG::LIB::UnordSet& a_v,
    SBG::LIB::UnordSet& b_v);


// auto return type we’ll let the compiler deduce what the return type should be from the return statement
template<typename M>
auto max_diff(M& cost_matrix, SBG::LIB::UnordSet& partition_a, SBG::LIB::UnordSet& partition_b, const WeightedSBGraph& graph)
{
    // cost_matrix is sort by gain, so the first is the maximum gain
    auto g = cost_matrix.begin();

    auto gain_object = *g;
#if KERNIGHAN_LIN_SBG_DEBUG
    cout << "The best is " << *g << endl;
#endif

    // remove it, we need to update those values that
    cost_matrix.erase(g);

    return gain_object;
}


std::ostream& operator<<(std::ostream& os, const GainObject& gain);

std::ostream& operator<<(std::ostream& os, const CostMatrix& cost_matrix);

std::ostream& operator<<(std::ostream& os, const KLBipartResult& result);

std::ostream& operator<<(std::ostream& os, const kl_sbg_partitioner_result& result);

}; // namespace sbg_partitioner
