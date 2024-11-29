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

#include <algorithm>
#include <chrono>
#include <future>
#include <map>
#include <rapidjson/document.h>
#include <rapidjson/filewritestream.h>
#include <rapidjson/prettywriter.h>
#include <rapidjson/writer.h>
#include <set>
#include <util/logger.hpp>

#include "build_sb_graph.hpp"
#include "kernighan_lin_partitioner.hpp"


#define PARTITION_IMBALANCE_DEBUG 1
#define PARTITION_IMBALANCE_PROFILE 1


// This code is based on https://github.com/CIFASIS/sbg-partitioner/discussions/17


using namespace std;

using namespace SBG::LIB;
using namespace SBG::Util;

namespace sbg_partitioner {

// Using unnamed namespace to define functions with internal linkage
namespace {

constexpr bool multithreading_enabled = false;

struct GainObjectImbalance {
    size_t i;
    size_t j;
    int gain;
    OrdSet ec_nodes_i;
    OrdSet ic_nodes_i;
    size_t size_i;
    OrdSet ec_nodes_j;
    OrdSet ic_nodes_j;
    size_t size_j;
};


struct KLBipartResult {
    SBG::LIB::OrdSet A;
    SBG::LIB::OrdSet B;
    int gain;
};


struct kl_sbg_partitioner_result
{
    size_t i;
    size_t j;
    int gain;
    SBG::LIB::OrdSet A;
    SBG::LIB::OrdSet B;
};


ostream& operator<<(ostream& os, const kl_sbg_partitioner_result& result)
{
    os << "{ partition results: "
       << result.i
       << ", "
       << result.j
       << ", "
       << result.gain
       << ", A: "
       << result.A
       << ", B: "
       << result.B
       << "}";

    return os;
}


template<typename G>
struct GainObjectComparatorTemplate {
    bool operator()(const G& gain_1, const G& gain_2) const
    {
        return gain_1.gain >= gain_2.gain;
    }
};


using ec_ic = std::pair<SBG::LIB::OrdPWMDInter , SBG::LIB::OrdPWMDInter>;

using GainObjectImbalanceComparator = GainObjectComparatorTemplate<GainObjectImbalance>;

using CostMatrixImbalance = std::set<GainObjectImbalance, GainObjectImbalanceComparator>;


ostream& operator<<(ostream& os, const KLBipartResult& result)
{
    os << "{ gain: " << result.gain << ", A: " << result.A << ", B: " << result.B << "}";

    return os;
}


ostream& operator<<(ostream& os, const GainObjectImbalance& gain)
{
    os << "< Node: ("
       << gain.i
       << ", size: "
       << gain.size_i
       << "), Node: ("
       << gain.j
       << ", size: "
       << gain.size_j
       << "), gain: "
       << gain.gain
       << " >";

    return os;
}


ostream& operator<<(ostream& os, const CostMatrixImbalance& cost_matrix)
{
    os << "{ ";
    for (const auto& o : cost_matrix) {
        os << o << " ";
    }
    os << " }";

    return os;
}


pair<unsigned, unsigned>
compute_lmin_lmax(const WeightedSBGraph& graph, unsigned number_of_partitions, const float imbalance_epsilon)
{
    unsigned w_v = get_node_size(graph.V(), graph.get_node_weights());
    unsigned B = ceil(w_v / number_of_partitions);
    int im = imbalance_epsilon * B;
    unsigned LMin = B - im;
    unsigned LMax = B + im;

    return make_pair(LMin, LMax);
}


size_t get_c_ab(
    const OrdSet& a, const OrdSet& b,
    const CanonPWMap& map_1,
    const CanonPWMap& map_2,
    const EdgeCost& costs)
{
    auto f = [](auto& a, auto& b, const CanonPWMap& map_1, const CanonPWMap& map_2) {
        OrdSet comm_edges;
        auto d = preImage(a, map_1);
        auto im = image(d, map_2);
        auto inters = intersection(b, im);
        auto edges = preImage(inters, map_2);
        edges = intersection(edges, d);
        comm_edges = cup(edges, comm_edges);

        return comm_edges;
    };

    auto intersection1 = f(a, b, map_1, map_2);
    auto intersection2 = f(a, b, map_2, map_1);

    auto communication_edges = cup(intersection1, intersection2);

    cout << "Comm between " << a << " and " << b << ": " << communication_edges << endl;

    size_t comm_size = get_edge_set_cost(communication_edges, costs);

    return comm_size;
}


#if PARTITION_IMBALANCE_PROFILE
static chrono::milliseconds total_duration;
#endif
ec_ic compute_EC_IC(
    const OrdSet& partition,
    const OrdSet& nodes,
    const OrdSet& partition_2,
    const CanonPWMap& map_1,
    const CanonPWMap& map_2)
{
    OrdSet ec, ic;
    auto d = preImage(nodes, map_1);
    auto im = image(d, map_2);
    auto ic_nodes = intersection(partition, im);
    ic_nodes = difference(ic_nodes, nodes);
    auto ec_nodes = difference(im, ic_nodes);
    ec_nodes = intersection(ec_nodes, partition_2);
    auto ic_ = intersection(preImage(ic_nodes, map_2), d);
    auto ec_ = intersection(preImage(ec_nodes, map_2), d);
    ec = cup(ec_, ec);
    ic = cup(ic_, ic);

    return make_pair(ec, ic);
}


unsigned get_imbalance_size(unsigned min_imbal_part, unsigned max_imbal_part, unsigned size_node_a, unsigned size_node_b, unsigned current_moved_size)
{
    cout << "nodes_imbal_part = " << max_imbal_part << ", " << size_node_a << endl;
    unsigned nodes_imbal_part = size_node_a > 0 ? floor(max_imbal_part / size_node_a) : 0;

    unsigned min_nodes_imbal_part = size_node_b > 0 ? floor(min_imbal_part / size_node_b) : 0;

    if (nodes_imbal_part > min_imbal_part) {
        nodes_imbal_part = min_nodes_imbal_part;
    }

    unsigned remaining_node_a = unsigned(size_node_a - current_moved_size);
    cout << "min between " << remaining_node_a << " and " << current_moved_size << endl;
    nodes_imbal_part = std::min(unsigned(current_moved_size), remaining_node_a);
    unsigned new_size_a = current_moved_size + nodes_imbal_part;

    return new_size_a;
}


GainObjectImbalance get_gain(
    unsigned idx_a,
    OrdSet& nodes_a,
    const OrdSet& partition_a,
    unsigned size_a,
    unsigned idx_b,
    OrdSet& nodes_b,
    const OrdSet& partition_b,
    unsigned size_b,
    const WeightedSBGraph& graph,
    const NodeWeight& node_weight)
{
    OrdSet a, b, rest_a, rest_b;
    tie(a, rest_a) = cut_interval_by_dimension(nodes_a, node_weight, size_a);
    tie(b, rest_b) = cut_interval_by_dimension(nodes_b, node_weight, size_b);

    // Now, compute external and internal cost for both maps
    OrdSet ec_nodes_a_1, ic_nodes_a_1;
    tie(ec_nodes_a_1, ic_nodes_a_1) = compute_EC_IC(partition_a, a, partition_b, graph.map1(), graph.map2());

    OrdSet ec_nodes_a_2, ic_nodes_a_2;
    tie(ec_nodes_a_2, ic_nodes_a_2) = compute_EC_IC(partition_a, a, partition_b, graph.map2(), graph.map1());

    // Get the union between both external and internal costs for both combination of maps
    OrdSet ec_nodes_a, ic_nodes_a;
    ec_nodes_a = cup(ec_nodes_a_1, ec_nodes_a_2);
    ic_nodes_a = cup(ic_nodes_a_1, ic_nodes_a_2);

    cout << "Node " << idx_a << ", " << nodes_a << " ec: " << ec_nodes_a << " and ic: " << ic_nodes_a << endl;

    size_t ec_a = get_edge_set_cost(ec_nodes_a, graph.get_edge_costs());
    size_t ic_a = get_edge_set_cost(ic_nodes_a, graph.get_edge_costs());
    int d_a = ec_a - ic_a;

    // Same as before for partition b
    OrdSet ec_nodes_b_1, ic_nodes_b_1;
    tie(ec_nodes_b_1, ic_nodes_b_1) = compute_EC_IC(partition_b, b, partition_a, graph.map1(), graph.map2());

    OrdSet ec_nodes_b_2, ic_nodes_b_2;
    tie(ec_nodes_b_2, ic_nodes_b_2) = compute_EC_IC(partition_b, b, partition_a, graph.map2(), graph.map1());

    OrdSet ec_nodes_b, ic_nodes_b;
    ec_nodes_b = cup(ec_nodes_b_1, ec_nodes_b_2);
    ic_nodes_b = cup(ic_nodes_b_1, ic_nodes_b_2);

    // cout << "Node: " << idx_b << ", " << nodes_b << " ec: " << ec_nodes_b << " and ic: " << ic_nodes_b << endl;

    size_t ec_b = get_edge_set_cost(ec_nodes_b, graph.get_edge_costs());
    size_t ic_b = get_edge_set_cost(ic_nodes_b, graph.get_edge_costs());
    int d_b = ec_b - ic_b;

    // Get communication between a and b
    size_t c_ab = get_c_ab(a, b, graph.map1(), graph.map2(), graph.get_edge_costs());

    // calculate gain
    int gain = d_a + d_b - 2 * c_ab;

    auto gain_obj = GainObjectImbalance{idx_a, idx_b, gain, ec_nodes_a, ic_nodes_a, size_a, ec_nodes_b, ic_nodes_b, size_b};

    return gain_obj;
}


void compute_exchange(unsigned i, unsigned j, OrdSet& partition_a, unsigned current_size_a,
    OrdSet& partition_b, unsigned current_size_b, const WeightedSBGraph& graph, const NodeWeight& node_weight,
    unsigned LMin, unsigned LMax, CostMatrixImbalance& cost_matrix)
{
    auto nodes_a = OrdSet(partition_a[i]);
    auto nodes_b = OrdSet(partition_b[j]);

    // Take the minimum partition size. We add 1 because it includes the last element
    size_t size_node_a = get_node_size(nodes_a, node_weight);
    size_t size_node_b = get_node_size(nodes_b, node_weight);
    size_t min_size = min(size_node_a, size_node_b);

    // No problem here, a is just a copy of partition_a[i], same for b
    GainObjectImbalance gain_obj = get_gain(i, nodes_a, partition_a, min_size, j, nodes_b, partition_b, min_size, graph, node_weight);

    // gain is greater than 0 and we are not moving all elements of node_a
    bool is_imbalance_enabled = LMin > 0 or LMax > 0;
    bool is_gain_positive = gain_obj.gain > 0;
    if (is_imbalance_enabled and is_gain_positive and size_node_a > min_size) {
        unsigned max_imbal_part = LMax > current_size_b ? LMax - current_size_b : 0; // this is what we can move to b

        unsigned min_imbal_part = current_size_a > LMin ? current_size_a - LMin : 0; // this is what we can move from a

        unsigned new_size_a = get_imbalance_size(min_imbal_part, max_imbal_part, size_node_a, size_node_b, min_size);

        GainObjectImbalance gain_obj_imbalance = get_gain(i, nodes_a, partition_a, new_size_a, j, nodes_b, partition_b, min_size, graph, node_weight);

        cout << "is gain better? " << gain_obj << ", " << gain_obj_imbalance << endl;

        if (gain_obj_imbalance.gain > gain_obj.gain) {
            gain_obj = move(gain_obj_imbalance);
        }
    }

    // gain is greater than 0 and we are not moving all elements of node_b
    if (is_imbalance_enabled and is_gain_positive and size_node_b > min_size) {
        unsigned max_imbal_part = LMax > current_size_a ? LMax - current_size_a : 0; // this is what we can move to a

        unsigned min_imbal_part = current_size_b > LMin ? current_size_b - LMin : 0; // this is what we can move from b

        unsigned new_size_b = get_imbalance_size(min_imbal_part, max_imbal_part, size_node_b, size_node_a, min_size);

        GainObjectImbalance gain_obj_imbalance = get_gain(i, nodes_a, partition_a, min_size, j, nodes_b, partition_b, new_size_b, graph, node_weight);

        cout << "is gain better? " << gain_obj << ", " << gain_obj_imbalance << endl;

        if (gain_obj_imbalance.gain > gain_obj.gain) {
            gain_obj = move(gain_obj_imbalance);
        }
    }

    cost_matrix.emplace(move(gain_obj));
}


CostMatrixImbalance generate_gain_matrix(
    const WeightedSBGraph& graph,
    const NodeWeight& node_weight,
    OrdSet& partition_a,
    OrdSet& partition_b,
    unsigned LMin,
    unsigned LMax)
{
    CostMatrixImbalance cost_matrix;

    unsigned p_size_a = get_node_size(partition_a, graph.get_node_weights());
    unsigned p_size_b = get_node_size(partition_b, graph.get_node_weights());

    for (size_t i = 0; i < partition_a.pieces().size(); i++) {
        for (size_t j = 0; j < partition_b.pieces().size(); j++) {
            compute_exchange(i, j, partition_a, p_size_a, partition_b, p_size_b, graph, node_weight, LMin, LMax, cost_matrix);
        }
    }

    return cost_matrix;
}


// Partition a and b (A_c and B_c in the definition) are the remining nodes to be visited, not the actual partitions
pair<pair<OrdSet, OrdSet>, pair<OrdSet, OrdSet>> update_sets(
    OrdSet& partition_a,
    OrdSet& partition_b,
    OrdSet& current_moved_partition_a,
    OrdSet& current_moved_partition_b,
    const GainObjectImbalance& gain_object,
    const WeightedSBGraph& graph)
{
    auto node_a = OrdSet(partition_a[gain_object.i]);
    size_t partition_size_a = get_node_size(node_a, graph.get_node_weights());
    bool node_a_is_fully_used = partition_size_a == gain_object.size_i;
    OrdSet rest_a;
    if (not node_a_is_fully_used) {
        tie(node_a, rest_a) = cut_interval_by_dimension(node_a, graph.get_node_weights(), gain_object.size_i);
        cout << "cut_interval_by_dimension " << gain_object.size_i << ": " << node_a << rest_a << endl;
    }

    auto node_b = OrdSet(partition_b[gain_object.j]);
    size_t partition_size_b = get_node_size(node_b, graph.get_node_weights());
    bool node_b_is_fully_used = partition_size_b == gain_object.size_j;
    OrdSet rest_b;
    if (not node_b_is_fully_used) {
        tie(node_b, rest_b) = cut_interval_by_dimension(node_b, graph.get_node_weights(), gain_object.size_j);
        cout << "cut_interval_by_dimension " << gain_object.size_j << ": " << node_b << rest_b << endl;
    }
    cout << "we remove " << node_a << " from " << partition_a << " and we get: ";
    partition_a = difference(partition_a, node_a);
    cout << partition_a << endl;
    cout << "we remove " << node_b << " from " << partition_b << " and we get: ";
    partition_b = difference(partition_b, node_b);
    cout << partition_b << endl;

    current_moved_partition_a = cup(current_moved_partition_a, node_a);
    current_moved_partition_b = cup(current_moved_partition_b, node_b);

    return make_pair(make_pair(node_a, rest_a), make_pair(node_b, rest_b));
}


void update_diff(
    CostMatrixImbalance& cost_matrix,
    OrdSet& remaining_partition_a,
    OrdSet& moved_from_partition_a,
    pair<OrdSet, OrdSet> affected_node_a,
    OrdSet& remaining_partition_b,
    OrdSet& moved_from_partition_b,
    pair<OrdSet, OrdSet> affected_node_b,
    const WeightedSBGraph& graph,
    const NodeWeight& node_weight,
    const GainObjectImbalance& gain_object,
    unsigned LMin,
    unsigned LMax)
{
    cout << affected_node_a.first << ", " << affected_node_a.second << endl;
    cout << affected_node_b.first << ", " << affected_node_b.second << endl;

    // Firstly, check if indexes need fixing. Three possible causes.
    size_t affected_node_a_size = get_node_size(affected_node_a.second, node_weight);
    bool node_a_fully_used = affected_node_a_size == 0;

    size_t affected_node_b_size = get_node_size(affected_node_b.second, node_weight);
    bool node_b_fully_used = affected_node_b_size == 0;

    unsigned size_a = get_node_size(remaining_partition_a, node_weight);
    size_a += get_node_size(moved_from_partition_b, node_weight);

    unsigned size_b = get_node_size(remaining_partition_b, node_weight);
    size_b += get_node_size(moved_from_partition_a, node_weight);

    // all the interval was used
    // fix indexes:
    if (node_a_fully_used or node_b_fully_used) {
        CostMatrixImbalance new_cost_matrix;
        for (auto g : cost_matrix) {
            if (node_a_fully_used and g.i == gain_object.i) {
                continue;
            }

            if (node_b_fully_used and g.j == gain_object.j) {
                continue;
            }

            if (node_a_fully_used and g.i > gain_object.i) {
                g.i--;
            }

            if (node_b_fully_used and g.j > gain_object.j) {
                g.j--;
            }

            new_cost_matrix.insert(g);
        }
        cost_matrix = new_cost_matrix;
    }

    if (not node_a_fully_used) {
        CostMatrixImbalance new_cost_matrix;
        for (auto g : cost_matrix) {
            if (g.i == gain_object.i) {
                compute_exchange(gain_object.i, g.j, remaining_partition_a, size_a, remaining_partition_b, size_b, graph, node_weight, LMin, LMax, new_cost_matrix);
            } else {
                new_cost_matrix.insert(g);
            }
        }

        cost_matrix = new_cost_matrix;
    }

    if (not node_b_fully_used) {
        CostMatrixImbalance new_cost_matrix;
        for (auto g : cost_matrix) {
            if (g.j == gain_object.j) {
                compute_exchange(g.i, gain_object.j, remaining_partition_a, size_a, remaining_partition_b, size_b, graph, node_weight, LMin, LMax, new_cost_matrix);
            } else {
                new_cost_matrix.insert(g);
            }
        }

        cost_matrix = new_cost_matrix;
    }

    CostMatrixImbalance new_cost_matrix;
    for (auto g : cost_matrix) {
        bool change = false;
        if (not isEmpty(intersection(g.ic_nodes_i, gain_object.ic_nodes_i)) or not isEmpty(intersection(g.ec_nodes_i, gain_object.ec_nodes_j))) {
            g.ic_nodes_i = difference(g.ic_nodes_i, gain_object.ic_nodes_i);
            g.ec_nodes_i = difference(g.ec_nodes_i, gain_object.ec_nodes_j);
            change = true;
        }

        if (not isEmpty(intersection(g.ic_nodes_j, gain_object.ic_nodes_j)) or not isEmpty(intersection(g.ec_nodes_j, gain_object.ec_nodes_i))) {
            g.ic_nodes_j = difference(g.ic_nodes_j, gain_object.ic_nodes_j);
            g.ec_nodes_j = difference(g.ec_nodes_j, gain_object.ec_nodes_i);
            change = true;
        }

        if (change) {
            size_t ec_i = get_edge_set_cost(g.ec_nodes_i, graph.get_edge_costs());
            size_t ic_i = get_edge_set_cost(g.ic_nodes_i, graph.get_edge_costs());
            int d_i = ec_i - ic_i;

            size_t ec_j = get_edge_set_cost(g.ec_nodes_j, graph.get_edge_costs());
            size_t ic_j = get_edge_set_cost(g.ic_nodes_j, graph.get_edge_costs());
            int d_j = ec_j - ic_j;

            // Get communication between a and b
            size_t c_ab = get_c_ab(remaining_partition_a[g.i], remaining_partition_b[g.j], graph.map1(), graph.map2(), graph.get_edge_costs());

            // calculate gain
            int gain = d_i + d_j - 2 * c_ab;
            g.gain = gain;
        }


        new_cost_matrix.insert(g);
    }
    cost_matrix = new_cost_matrix;

#if PARTITION_IMBALANCE_DEBUG
    cout << remaining_partition_a << ", " << remaining_partition_b << ", " << gain_object << ", " << cost_matrix << endl;
#endif
}


// auto return type we’ll let the compiler deduce what the return type should be from the return statement
template<typename M>
auto max_diff(M& cost_matrix, SBG::LIB::OrdSet& partition_a, SBG::LIB::OrdSet& partition_b, const WeightedSBGraph& graph)
{
    // cost_matrix is sort by gain, so the first is the maximum gain
    auto g = cost_matrix.begin();

    auto gain_object = *g;
#if PARTITION_IMBALANCE_DEBUG
    cout << "The best is " << *g << endl;
#endif

    // remove it, we need to update those values that
    cost_matrix.erase(g);

    return gain_object;
}


void update_sum(
    int& par_sum,
    int g,
    int& max_par_sum,
    pair<OrdSet, OrdSet>& max_par_sum_set,
    const OrdSet& a_v,
    const OrdSet& b_v)
{
    par_sum += g;
    if (par_sum > max_par_sum) {
        max_par_sum = par_sum;
        max_par_sum_set = make_pair(a_v, b_v);
    }
}


int kl_sbg_imbalance(
    const WeightedSBGraph& graph,
    OrdSet& partition_a,
    OrdSet& partition_b,
    unsigned LMin,
    unsigned LMax)
{
#if PARTITION_IMBALANCE_DEBUG
    cout << "Algorithm starts with " << partition_a << ", " << partition_b << endl;
#endif
    auto a_c = partition_a;
    auto b_c = partition_b;
    int max_par_sum = 0;
    auto max_par_sum_set = make_pair(OrdSet(), OrdSet());
    int par_sum = 0;
    OrdSet a_v = OrdSet();
    OrdSet b_v = OrdSet();
    const auto node_weights = graph.get_node_weights();

    CostMatrixImbalance gm = generate_gain_matrix(graph, node_weights, partition_a, partition_b, LMin, LMax);

#if PARTITION_IMBALANCE_DEBUG
        cout << LMin << ", "
             << LMax
             << gm << endl;
#endif

    while ((not isEmpty(a_c)) and (not isEmpty(b_c))) {
        cout << "inside the while " << a_c << b_c << endl;
        GainObjectImbalance g = max_diff(gm, a_c, b_c, graph);
        cout << g << endl;
        pair<OrdSet, OrdSet> a_, b_;
        tie(a_, b_) = update_sets(a_c, b_c, a_v, b_v, g, graph);
        update_diff(gm, a_c, a_v, a_, b_c, b_v, b_, graph, node_weights, g, LMin, LMax);
        update_sum(par_sum, g.gain, max_par_sum, max_par_sum_set, a_v, b_v);
    }

    if (max_par_sum > 0) {

        partition_a = cup(difference(partition_a, max_par_sum_set.first), max_par_sum_set.second);
        partition_b = cup(difference(partition_b, max_par_sum_set.second), max_par_sum_set.first);

        flatten_set(partition_a, graph);
        flatten_set(partition_b, graph);
    }

#if PARTITION_IMBALANCE_DEBUG
    cout << "so it ends with " << max_par_sum << ", " << partition_a << ", " << partition_b << endl;
#endif
    return max_par_sum;
}


KLBipartResult kl_sbg_bipart_imbalance(const WeightedSBGraph& graph, OrdSet& partition_a,
    OrdSet& partition_b, unsigned LMin, unsigned LMax)
{
    int gain = kl_sbg_imbalance(graph, partition_a, partition_b, LMin, LMax);

#if PARTITION_IMBALANCE_DEBUG
    cout << "Final: " << partition_a << ", " << partition_b << endl;
#endif

    return KLBipartResult{partition_a, partition_b, gain};
}


kl_sbg_partitioner_result kl_sbg_partitioner_function(
    const WeightedSBGraph& graph, PartitionMap& partitions, unsigned LMin, unsigned LMax,
    vector<kl_sbg_partitioner_result>& gains)
{
    map<size_t, OrdSet> adjacents;
    kl_sbg_partitioner_result best_gain = kl_sbg_partitioner_result{ 0, 0, -1, OrdSet(), OrdSet()};
    for (size_t i = 0; i < partitions.size(); i++) {

        if (adjacents.find(i) == adjacents.end()) {
            adjacents[i] = get_adjacents(graph, partitions[i]);
        }

        for (size_t j = i + 1; j < partitions.size(); j++) {

            if (isEmpty(intersection(adjacents[i], partitions[j]))) {
                cout << "No connections between " << partitions[i] << " and " << partitions[j] << " is empty" << endl;
                continue;
            }

            auto gain_comp = [i, j] (const kl_sbg_partitioner_result& g) {
                return (g.i == i and g.j == j) or (g.i == j and g.j == i);
            };

            auto gain_it = find_if(gains.begin(), gains.end(), gain_comp);
            if (gain_it != gains.end()) {
                cout << "Between " << i << " and " << j << " was already computed, " << *gain_it  << endl;
                continue;
            }

            auto p_1_copy = partitions[i];
            auto p_2_copy = partitions[j];
            KLBipartResult current_gain = kl_sbg_bipart_imbalance(graph, p_1_copy, p_2_copy, LMin, LMax);
    #if PARTITION_IMBALANCE_DEBUG
            cout << "current_gain " << current_gain << endl;
    #endif
            gains.emplace_back(kl_sbg_partitioner_result{ i, j, current_gain.gain, current_gain.A, current_gain.B });
        }
    }

    for_each(gains.begin(), gains.end(), [&best_gain](const kl_sbg_partitioner_result& current_gain) {
        if (current_gain.gain > best_gain.gain) {
            best_gain = current_gain;
        }
    });

    return best_gain;
}


kl_sbg_partitioner_result kl_sbg_partitioner_multithreading(
    const WeightedSBGraph& graph, PartitionMap& partitions, unsigned LMin, unsigned LMax,
    vector<kl_sbg_partitioner_result>& gains)
{
    kl_sbg_partitioner_result best_gain = kl_sbg_partitioner_result{ 0, 0, -1, OrdSet(), OrdSet()};
    vector<future<kl_sbg_partitioner_result>> workers;
    map<size_t, OrdSet> adjacents;
    for (size_t i = 0; i < partitions.size(); i++) {

        if (adjacents.find(i) == adjacents.end()) {
            adjacents[i] = get_adjacents(graph, partitions[i]);
        }

        for (size_t j = i + 1; j < partitions.size(); j++) {

            if (isEmpty(intersection(adjacents[i], partitions[j]))) {
                cout << "No connections between " << partitions[i] << " and " << partitions[j] << " is empty" << endl;
                continue;
            }

            auto gain_comp = [i, j] (const kl_sbg_partitioner_result& g) {
                return (g.i == i and g.j == j) or (g.i == j and g.j == i);
            };

            auto gain_it = find_if(gains.begin(), gains.end(), gain_comp);
            if (gain_it != gains.end()) {
                cout << "Between " << i << " and " << j << " was already computed, " << *gain_it  << endl;
                continue;
            }

            auto th = async([&graph, &partitions, i, j, LMin, LMax] () {
                OrdSet p_1_copy = partitions[i];
                OrdSet p_2_copy = partitions[j];
                KLBipartResult results = kl_sbg_bipart_imbalance(graph, p_1_copy, p_2_copy, LMin, LMax);
                return kl_sbg_partitioner_result{i, j, results.gain, results.A, results.B};
            });
            workers.push_back(move(th));
        }
    }

    for_each(workers.begin(), workers.end(), [&best_gain, &gains] (future<kl_sbg_partitioner_result>& th) {
        // here we wait for each thread to finish and get its results
        auto current_gain = th.get();
        gains.emplace_back(current_gain);
    });

    for_each(gains.begin(), gains.end(), [&best_gain] (const kl_sbg_partitioner_result& current_gain) {
        if (current_gain.gain > best_gain.gain) {
            best_gain = current_gain;
        }
    });

    return best_gain;
}


void kl_sbg_imbalance_partitioner(
    const WeightedSBGraph& graph, PartitionMap& partitions, const float imbalance_epsilon)
{
    auto [LMin, LMax] = imbalance_epsilon > 0.0 ? compute_lmin_lmax(graph, partitions.size(), imbalance_epsilon) : make_pair<unsigned, unsigned>(0, 0);
    bool change = true;
    int counter = 0;

    vector<kl_sbg_partitioner_result> gains;
    while (change) {
        cout << "*****ITERATION NUMBER " << counter++ << endl;
        change = false;

        kl_sbg_partitioner_result best_gain;
        if (multithreading_enabled) {
            best_gain = kl_sbg_partitioner_multithreading(graph, partitions, LMin, LMax, gains);
        } else {
            best_gain = kl_sbg_partitioner_function(graph, partitions, LMin, LMax, gains);
        }

        cout << "Best gain results is: " << best_gain << endl;

        auto gain_comp = [&best_gain](const kl_sbg_partitioner_result& g) {
            return g.i == best_gain.i or g.j == best_gain.j
                or g.i == best_gain.j or g.j == best_gain.i;
        };

        // now, apply changes
        constexpr int strategy = 2;
        // first strategy
        switch (strategy)
        {
        case 1:
            if (best_gain.gain > 0) {
                change = true;
                partitions[best_gain.i] = best_gain.A;
                partitions[best_gain.j] = best_gain.B;

                gains.erase(std::remove_if(gains.begin(), gains.end(), gain_comp), gains.end());
            }

            break;
        default:

            int it_counter = 0;
            while (not gains.empty() and best_gain.gain > 0) {
                cout << "change number " << it_counter << " changing " << best_gain.i << ", " << best_gain.j << endl;
                it_counter++;
                change = true;
                partitions[best_gain.i] = best_gain.A;
                partitions[best_gain.j] = best_gain.B;

                gains.erase(std::remove_if(gains.begin(), gains.end(), gain_comp), gains.end());

                cout << "best gain is " << best_gain << endl;
                cout << "and vector is ";
                for_each(gains.begin(), gains.end(), [](const kl_sbg_partitioner_result& g) { cout << g << " "; });
                cout << endl;

                if (not gains.empty()){
                    auto max_gain_it = max_element(gains.begin(), gains.end(),
                        [] (const kl_sbg_partitioner_result& a, const kl_sbg_partitioner_result& b) {
                            return a.gain < b.gain;
                        });
                    if (max_gain_it != gains.end()) {
                        best_gain = *max_gain_it;
                    } else {
                        best_gain = kl_sbg_partitioner_result{ 0, 0, -Inf, {}, {}} ;
                    }
                }
            }

            break;
        }

        for (size_t i = 0; i < partitions.size(); i++) {
            SBG_LOG << i << ": " << partitions[i] << endl;
        }
    }

    for (size_t i = 0; i < partitions.size(); i++) {
        SBG_LOG << i << ": " << partitions[i] << endl;
    }
}


string get_pretty_sb_graph(const SBG::LIB::CanonSBG& g)
{
    rapidjson::Document json_doc;
    rapidjson::Document::AllocatorType& allocator = json_doc.GetAllocator();
    json_doc.SetObject();

    rapidjson::Value obj_nodes(rapidjson::kArrayType);
    for (size_t i = 0; i < g.V().size(); i++) {
        const SetPiece& set_piece = g.V()[i];
        rapidjson::Value obj_set_piece_array(rapidjson::kArrayType);
        for (const Interval& interval : set_piece.intervals()) {
            rapidjson::Value obj_interval_array(rapidjson::kArrayType);

            rapidjson::Value begin(rapidjson::kNumberType);
            begin.SetUint(interval.begin());
            obj_interval_array.PushBack(begin, allocator);

            rapidjson::Value end(rapidjson::kNumberType);
            end.SetUint(interval.end());
            obj_interval_array.PushBack(end, allocator);

            obj_set_piece_array.PushBack(obj_interval_array, allocator);
        }

        obj_nodes.PushBack(obj_set_piece_array, allocator);
    }

    json_doc.AddMember("nodes", obj_nodes, allocator);

    auto map_parser = [&json_doc, &allocator] (const auto& maps, const char* key) {
        rapidjson::Value obj_maps(rapidjson::kArrayType);
        for (const auto& map1 : maps) {
            rapidjson::Value map_1_obj(rapidjson::kObjectType);

            const auto& domain = map1.dom()[0];
            rapidjson::Value domain_obj(rapidjson::kArrayType);
            for (const auto& interval : domain.intervals()) {
                rapidjson::Value obj_interval_array(rapidjson::kArrayType);

                rapidjson::Value begin(rapidjson::kNumberType);
                begin.SetUint(interval.begin());
                obj_interval_array.PushBack(begin, allocator);

                rapidjson::Value end(rapidjson::kNumberType);
                end.SetUint(interval.end());
                obj_interval_array.PushBack(end, allocator);

                domain_obj.PushBack(obj_interval_array, allocator);
            }
            map_1_obj.AddMember("domain", domain_obj, allocator);

            rapidjson::Value exp_obj(rapidjson::kArrayType);
            for (const auto& exp : map1.exp()) {
                rapidjson::Value obj_exp_array(rapidjson::kArrayType);

                rapidjson::Value slope(rapidjson::kNumberType);
                float _slope = exp.slope().numerator();
                slope.SetFloat(_slope);
                obj_exp_array.PushBack(slope, allocator);

                rapidjson::Value offset(rapidjson::kNumberType);
                float _offset = exp.offset().numerator();
                slope.SetFloat(_offset);
                obj_exp_array.PushBack(offset, allocator);

                exp_obj.PushBack(obj_exp_array, allocator);
            }
            map_1_obj.AddMember("exp", exp_obj, allocator);

            obj_maps.PushBack(map_1_obj, allocator);
        }

        rapidjson::Value key_value(key, allocator);

        json_doc.AddMember(key_value, obj_maps, allocator);
    };

    map_parser(g.map1().maps(), "map1");
    map_parser(g.map2().maps(), "map2");

    // Write the JSON data to the file
    rapidjson::StringBuffer s;
    rapidjson::Writer<rapidjson::StringBuffer> writer(s);
    json_doc.Accept(writer);

    string json_data = string(s.GetString());

    return json_data;
}


}


string partitionate_nodes(
    const std::string& filename,
    const unsigned number_of_partitions,
    const float epsilon,
    optional<string>& graph_str)
{
    auto sb_graph = build_sb_graph(filename.c_str());

    cout << sb_graph << endl;
    cout << "sb graph created!" << endl;

    auto partitions = best_initial_partition(sb_graph, number_of_partitions);

    kl_sbg_imbalance_partitioner(sb_graph, partitions, epsilon);

    // // just for debugging
    cout << endl;
    for (unsigned i = 0; i < partitions.size(); i++) {
        cout << i << " " << partitions[i] << endl;
    }
    cout << endl;

    sanity_check(sb_graph, partitions, number_of_partitions);

    if (graph_str){
        graph_str = get_pretty_sb_graph(sb_graph);
    }

    string output = get_output(partitions);

    return output;
}


pair<WeightedSBGraph, PartitionMap> partitionate_nodes_for_metrics(
    const string& filename,
    const unsigned number_of_partitions,
    const float epsilon)
{
    auto sb_graph = build_sb_graph(filename.c_str());
    // auto sb_graph = create_air_conditioners_graph();

    cout << sb_graph << endl;
    cout << "sb graph created!" << endl;

    auto partitions = best_initial_partition(sb_graph, number_of_partitions);

    kl_sbg_imbalance_partitioner(sb_graph, partitions, epsilon);

    return {sb_graph, partitions};
}

}