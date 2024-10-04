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

#include <future>
#include <set>
#include <sbg/sbg.hpp>
#include <thread>
#include <util/logger.hpp>
#include <utility>
#include <vector>

#include <unistd.h>

#include "build_sb_graph.hpp"
#include "kernighan_lin_sbg.hpp"


#define KERNIGHAN_LIN_SBG_DEBUG 1


// This code is based on https://github.com/CIFASIS/sbg-partitioner/discussions/17


using namespace std;

using namespace SBG::LIB;
using namespace SBG::Util;

namespace sbg_partitioner {


ostream& operator<<(ostream& os, const GainObject& gain)
{
    os << "< Node: ("
       << gain.i
       << ", "
       << gain.j
       << "), gain: "
       << gain.gain
       << ", size: "
       << gain.size
       << " >";

    return os;
}

std::ostream& operator<<(std::ostream& os, const CostMatrix& cost_matrix)
{
    os << "{ ";
    for (const auto& o : cost_matrix) {
        os << o << " ";
    }
    os << " }";

    return os;
}







// Consider exposing these functions to facilitate being reused



static GainObject compute_diff(
    size_t i, size_t j,
    UnordSet& partition_a,
    UnordSet& partition_b,
    const WeightedSBGraph& graph)
{
    // Firstly, copy both partitions
    auto a = UnordSet(partition_a[i]);
    auto b = UnordSet(partition_b[j]);

    // Take the minimum partition size. We add 1 because it includes the last element
    size_t size_node_a = get_node_size(a, graph.get_node_weights());
    size_t size_node_b = get_node_size(b, graph.get_node_weights());
    size_t min_size = min(size_node_a, size_node_b);

    int weight_a = get_set_cost(partition_a[i], graph.get_node_weights());
    int weight_b = get_set_cost(partition_b[j], graph.get_node_weights());

    if (min_size < size_t(weight_a) or min_size < size_t(weight_b)) {
        return GainObject{i, j, 0, 0};
    }

    // No problem here, a is just a copy of partition_a[i], same for b
    // We substract 1 because it includes the last element
    UnordSet rest_a, rest_b;
    tie(a, rest_a) = cut_interval_by_dimension(a, graph.get_node_weights(), min_size);
    tie(b, rest_b) = cut_interval_by_dimension(b, graph.get_node_weights(), min_size);

    // Now, compute external and internal cost for both maps
    UnordSet ec_nodes_a_1, ic_nodes_a_1;
    tie(ec_nodes_a_1, ic_nodes_a_1) = compute_EC_IC(partition_a, a, partition_b, graph.map1(), graph.map2());

    UnordSet ec_nodes_a_2, ic_nodes_a_2;
    tie(ec_nodes_a_2, ic_nodes_a_2) = compute_EC_IC(partition_a, a, partition_b, graph.map2(), graph.map1());

    // Get the union between both external and internal costs for both combination of maps
    UnordSet ec_nodes_a, ic_nodes_a;
    ec_nodes_a = cup(ec_nodes_a_1, ec_nodes_a_2);
    ic_nodes_a = cup(ic_nodes_a_1, ic_nodes_a_2);

    size_t ec_a = get_edge_set_cost(ec_nodes_a, graph.get_edge_costs());
    size_t ic_a = get_edge_set_cost(ic_nodes_a, graph.get_edge_costs());
    int d_a = ec_a - ic_a;

    // Same as before for partition b
    UnordSet ec_nodes_b_1, ic_nodes_b_1;
    tie(ec_nodes_b_1, ic_nodes_b_1) = compute_EC_IC(partition_b, b, partition_a, graph.map1(), graph.map2());

    UnordSet ec_nodes_b_2, ic_nodes_b_2;
    tie(ec_nodes_b_2, ic_nodes_b_2) = compute_EC_IC(partition_b, b, partition_a, graph.map2(), graph.map1());

    UnordSet ec_nodes_b, ic_nodes_b;
    ec_nodes_b = cup(ec_nodes_b_1, ec_nodes_b_2);
    ic_nodes_b = cup(ic_nodes_b_1, ic_nodes_b_2);

    size_t ec_b = get_edge_set_cost(ec_nodes_b, graph.get_edge_costs());
    size_t ic_b = get_edge_set_cost(ic_nodes_b, graph.get_edge_costs());
    int d_b = ec_b - ic_b;

    // Get communication between a and b
    size_t c_ab = get_c_ab(a, b, graph.map1(), graph.map2(), graph.get_edge_costs());

    // calculate gain
    int gain = d_a + d_b - 2 * c_ab;

    // Create Gain object
    const auto gain_obj = GainObject{i, j, gain, min_size};
    return gain_obj;
}


static CostMatrix generate_gain_matrix(
    UnordSet& partition_a,
    UnordSet& partition_b,
    const WeightedSBGraph& graph)
{
    CostMatrix cost_matrix;

    for (size_t i = 0; i < partition_a.pieces().size(); i++) {
        for (size_t j = 0; j < partition_b.pieces().size(); j++) {
            // Get gain object for both nodes
            auto gain_obj = compute_diff(i, j, partition_a, partition_b, graph);

            // Insert it to the cost matrix
            cost_matrix.insert(gain_obj);
        }
    }

    return cost_matrix;
}


static void remove_nodes_from_cost_matrix(
    CostMatrix& cost_matrix,
    UnordSet& partition,
    const GainObject& gain_object)
{
    vector<CostMatrix::iterator> to_be_removed;

    // delete each combination with partition i
    for (auto it = cost_matrix.begin(); it != cost_matrix.end(); ++it) {
        if (it->j == gain_object.j) {
            to_be_removed.push_back(it);
        }
    }

    for (auto& it : to_be_removed) {
        cost_matrix.erase(it);
    }
}


static void update_cost_matrix(
    UnordSet& partition_a,
    UnordSet& partition_b,
    const GainObject& gain_object,
    const WeightedSBGraph& graph,
    CostMatrix& cost_matrix)
{
    // Let's add gain objects to the cost matrix
    for (size_t i = 0; i < partition_a.size(); i++) {
        auto new_gain_object = compute_diff(i, gain_object.j, partition_a, partition_b, graph);
        cost_matrix.insert(new_gain_object);
    }
}


// Partition a and b (A_c and B_c in the definition) are the remining nodes to be visited, not the actual partitions
static pair<UnordSet, UnordSet> update_sets(
    UnordSet& partition_a,
    UnordSet& partition_b,
    UnordSet& current_moved_partition_a,
    UnordSet& current_moved_partition_b,
    const GainObject& gain_object,
    const NodeWeight& node_weight)
{
    auto node_a = UnordSet(partition_a[gain_object.i]);
    size_t partition_size_a = get_node_size(node_a, node_weight);
    bool node_a_is_fully_used = partition_size_a == gain_object.size;
    if (not node_a_is_fully_used) {
        UnordSet rest_a;
        tie(node_a, rest_a) = cut_interval_by_dimension(node_a, node_weight, gain_object.size);
    }

    auto node_b = UnordSet(partition_b[gain_object.j]);
    size_t partition_size_b = get_node_size(node_b, node_weight);
    bool node_b_is_fully_used = partition_size_b == gain_object.size;
    if (not node_b_is_fully_used) {
        UnordSet rest_b;
        tie(node_b, rest_b) = cut_interval_by_dimension(node_b, node_weight, gain_object.size);
    }

    // At least one of them should be fully used?
#if KERNIGHAN_LIN_SBG_DEBUG
    cout << boolalpha
         << "node_a_is_fully_used: "
         << node_a_is_fully_used
         << ", node_b_is_fully_used "
         << node_b_is_fully_used
         << endl;
#endif

    partition_a = difference(partition_a, node_a);
    partition_b = difference(partition_b, node_b);

    current_moved_partition_a = cup(current_moved_partition_a, node_a);
    current_moved_partition_b = cup(current_moved_partition_b, node_b);

    return make_pair(node_a, node_b);
}


static void update_diff(
    CostMatrix& cost_matrix,
    UnordSet& partition_a,
    UnordSet& partition_b,
    const WeightedSBGraph& graph,
    const GainObject& gain_object)
{
    if (false) { // we shoudl check on this
        remove_nodes_from_cost_matrix(cost_matrix, partition_a, gain_object);
        remove_nodes_from_cost_matrix(cost_matrix, partition_b, gain_object);

        update_cost_matrix(partition_b, partition_a, gain_object, graph, cost_matrix);
        update_cost_matrix(partition_a, partition_b, gain_object, graph, cost_matrix);
    }
    cost_matrix.clear();
    cost_matrix = generate_gain_matrix(partition_a, partition_b, graph);
#if KERNIGHAN_LIN_SBG_DEBUG
    cout << partition_a << ", " << partition_b << ", " << cost_matrix << endl;
#endif

}





int kl_sbg(const WeightedSBGraph& graph, UnordSet& partition_a, UnordSet& partition_b)
{
#if KERNIGHAN_LIN_SBG_DEBUG
    cout << "Algorithm starts with " << partition_a << ", " << partition_b << endl;
#endif
    auto a_c = partition_a;
    auto b_c = partition_b;
    int max_par_sum = 0;
    auto max_par_sum_set = make_pair(UnordSet(), UnordSet());
    int par_sum = 0;
    UnordSet a_v = UnordSet();
    UnordSet b_v = UnordSet();

    CostMatrix gm = generate_gain_matrix(partition_a, partition_b, graph);

#if KERNIGHAN_LIN_SBG_DEBUG
    cout << gm << endl;
#endif

    while ((not isEmpty(a_c)) and (not isEmpty(b_c))) {
        GainObject g = max_diff(gm, a_c, b_c, graph);
        UnordSet a_, b_;
        tie(a_, b_) = update_sets(a_c, b_c, a_v, b_v, g, graph.get_node_weights());
        update_diff(gm, a_c, b_c, graph, g);
        update_sum(par_sum, g.gain, max_par_sum, max_par_sum_set, a_v, b_v);
    }

    if (max_par_sum > 0) {
        partition_a = cup(difference(partition_a, max_par_sum_set.first), max_par_sum_set.second);
        partition_b = cup(difference(partition_b, max_par_sum_set.second), max_par_sum_set.first);
    }

#if KERNIGHAN_LIN_SBG_DEBUG
    cout << "so it ends with " << max_par_sum << ", " << partition_a << ", " << partition_b << endl;
#endif

    return max_par_sum;
}


KLBipartResult kl_sbg_bipart(const WeightedSBGraph& graph, SBG::LIB::UnordSet& partition_a, SBG::LIB::UnordSet& partition_b)
{
    auto partition_a_copy = partition_a;
    auto partition_b_copy = partition_b;
    int gain = kl_sbg(graph, partition_a_copy, partition_b_copy);

    while (not (partition_a_copy == partition_a) and not (partition_b_copy == partition_b)) {
        partition_a = partition_a_copy;
        partition_b = partition_b_copy;
        gain = max(kl_sbg(graph, partition_a_copy, partition_b_copy), gain);
#if KERNIGHAN_LIN_SBG_DEBUG
        cout << "gain: " << gain << endl;
#endif
    }

#if KERNIGHAN_LIN_SBG_DEBUG
    cout << "Final: " << partition_a << ", " << partition_b << endl;
#endif
    return KLBipartResult{partition_a, partition_b, gain};
}


ostream& operator<<(ostream& os, const kl_sbg_partitioner_result& result)
{
    os << "i: " << result.i << ": " << result.A << ", j: " << result.j << ", " << result.B << ", gain: " << result.gain;

    return os;
}


static kl_sbg_partitioner_result kl_sbg_partitioner_function(const WeightedSBGraph& graph, PartitionMap& partitions)
{
    kl_sbg_partitioner_result best_gain = kl_sbg_partitioner_result{ 0, 0, -1, UnordSet(), UnordSet()};
    for (size_t i = 0; i < partitions.size(); i++) {
        for (size_t j = i + 1; j < partitions.size(); j++) {
            auto p_1_copy = partitions[i];
            auto p_2_copy = partitions[j];
            KLBipartResult current_gain = kl_sbg_bipart(graph, p_1_copy, p_2_copy);
    #if KERNIGHAN_LIN_SBG_DEBUG
            cout << "current_gain " << current_gain << endl;
    #endif
            if (current_gain.gain > best_gain.gain) {
                best_gain.i = i;
                best_gain.j = j;
                best_gain.gain = current_gain.gain;
                best_gain.A = current_gain.A;
                best_gain.B = current_gain.B;
            }
        }
    }

    return best_gain;
}


static kl_sbg_partitioner_result kl_sbg_partitioner_multithreading(const WeightedSBGraph& graph, PartitionMap& partitions)
{
    kl_sbg_partitioner_result best_gain = kl_sbg_partitioner_result{ 0, 0, -1, UnordSet(), UnordSet()};
    vector<future<kl_sbg_partitioner_result>> workers;
    for (size_t i = 0; i < partitions.size(); i++) {
        for (size_t j = i + 1; j < partitions.size(); j++) {
            auto th = async([&graph, &partitions, i, j] () {
                UnordSet p_1_copy = partitions[i];
                UnordSet p_2_copy = partitions[j];
                KLBipartResult results = kl_sbg_bipart(graph, p_1_copy, p_2_copy);
                return kl_sbg_partitioner_result{i, j, results.gain, results.A, results.B};
            });
            workers.push_back(move(th));
        }
    }

    for_each(workers.begin(), workers.end(), [&best_gain] (future<kl_sbg_partitioner_result>& th) {
        // here we wait for each thread to finish and get its results
        auto current_gain = th.get();
        if (current_gain.gain > best_gain.gain) {
            best_gain = move(current_gain);
        }
    });

    return best_gain;
}


void kl_sbg_partitioner(const WeightedSBGraph& graph, PartitionMap& partitions)
{
    bool change = true;
    while (change) {
        change = false;

        kl_sbg_partitioner_result best_gain;
        if (multithreading_enabled) {
            best_gain = kl_sbg_partitioner_multithreading(graph, partitions);
        } else {
            best_gain = kl_sbg_partitioner_function(graph, partitions);
        }

        // now, apply changes
        if (best_gain.gain > 0) {
            change = true;
            partitions[best_gain.i] = best_gain.A;
            partitions[best_gain.j] = best_gain.B;
        }

        if (change and graph.V()[0].intervals().size() == 1) {
            flatten_set(partitions[best_gain.i], graph);
            flatten_set(partitions[best_gain.j], graph);
        }
    }

    for (size_t i = 0; i < partitions.size(); i++) {
        SBG_LOG << i << ": " << partitions[i] << endl;
    }
}

}
