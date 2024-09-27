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


#include "build_sb_graph.hpp"
#include "kernighan_lin_sbg.hpp"
#include "partition_imbalance.hpp"


#define PARTITION_IMBALANCE_DEBUG 1


using namespace std;

using namespace SBG::LIB;

namespace sbg_partitioner {


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


static void partition_imbalance(
    UnordSet& nodes,
    UnordSet& partition,
    UnordSet& partition_2,
    const NodeWeight& node_weights,
    vector<size_t>& i,
    vector<size_t>& i_2,
    unsigned LMin,
    unsigned LMax)
{
    unsigned p_size = get_node_size(partition, node_weights);
    unsigned max_imbal_part = LMax > p_size ? LMax - p_size : 0;
    cout << "For " << partition << " are " << LMax << ", " << p_size << ", " << max_imbal_part << endl;

    unsigned p_size_2 = get_node_size(partition_2, node_weights);
    unsigned min_imbal_part = p_size_2 > LMin ? p_size_2 - LMin : 0;
    cout << "For " << partition_2 << " are " << LMin << ", " << p_size_2 << ", " << min_imbal_part << endl;

    // assuming we only move one set piece
    assert(nodes.size() == 1);
    unsigned node_weight = get_set_cost(nodes[0], node_weights);
    cout << "node weight is " << node_weight << endl;

    unsigned nodes_imbal_part = floor(max_imbal_part / node_weight);

    unsigned min_nodes_imbal_part = floor(min_imbal_part / node_weight);

    if (nodes_imbal_part > min_imbal_part) {
        nodes_imbal_part = min_nodes_imbal_part;
    }

    if (nodes_imbal_part > 0 and unsigned(nodes_imbal_part) > get_node_size(nodes, NodeWeight())) {
        nodes_imbal_part = get_node_size(nodes, NodeWeight());
    }

    if (nodes_imbal_part > 0) {
        i.push_back(0);
        i_2.push_back(nodes_imbal_part);
    }
}


void compute_partition_imbalance(unsigned i, unsigned j, UnordSet& partition_a,
    UnordSet& partition_b, const WeightedSBGraph& graph, const NodeWeight& node_weight,
    int LMin, int LMax, CostMatrixImbalance& cost_matrix)
{
    vector<size_t> i_a;
    vector<size_t> i_b;

    auto nodes_a = UnordSet(partition_a[i]);
    auto nodes_b = UnordSet(partition_b[j]);

    // Take the minimum partition size. We add 1 because it includes the last element
    size_t size_node_a = get_node_size(nodes_a, node_weight);
    size_t size_node_b = get_node_size(nodes_b, node_weight);
    size_t min_size = min(size_node_a, size_node_b);

    // check how many nodes we can move from one side to another
    int weight_a = get_set_cost(partition_a[i], node_weight);
    int weight_b = get_set_cost(partition_b[j], node_weight);

    int nodes_bal_part_a = floor(min_size / weight_a);
    int nodes_bal_part_b = floor(min_size / weight_b);

    i_a.push_back(nodes_bal_part_a);
    i_b.push_back(nodes_bal_part_b);

    partition_imbalance(nodes_a, partition_a, partition_b, node_weight, i_a, i_b, LMin, LMax);

    partition_imbalance(nodes_b, partition_b, partition_a, node_weight, i_b, i_a, LMin, LMax);

    assert(i_a.size() == i_b.size());

    cout << partition_a << endl;
    cout << "i_a: ";
    for_each(i_a.cbegin(), i_a.cend(), [](const auto& s) { cout << s << " "; });
    cout << endl;

    cout << partition_b << endl;
    cout << "i_b: ";
    for_each(i_b.cbegin(), i_b.cend(), [](const auto& s) { cout << s << " "; });
    cout << endl;

    for (size_t idx = 0; idx < i_a.size(); idx++) {
        // No problem here, a is just a copy of partition_a[i], same for b
        // We substract 1 because it includes the last element
        UnordSet a, b, rest_a, rest_b;
        tie(a, rest_a) = cut_interval_by_dimension(nodes_a, node_weight, i_a[idx]);
        tie(b, rest_b) = cut_interval_by_dimension(nodes_b, node_weight, i_b[idx]);

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

        auto gain_obj = GainObjectImbalance{i, j, gain, i_a[idx], i_b[idx]};
        cost_matrix.emplace(move(gain_obj));
    }
}


CostMatrixImbalance generate_gain_matrix(
    const WeightedSBGraph& graph,
    const NodeWeight& node_weight,
    UnordSet& partition_a,
    UnordSet& partition_b,
    int LMin,
    int LMax)
{
    CostMatrixImbalance cost_matrix;

    for (size_t i = 0; i < partition_a.pieces().size(); i++) {
        for (size_t j = 0; j < partition_b.pieces().size(); j++) {
            compute_partition_imbalance(i, j, partition_a, partition_b, graph, node_weight, LMin, LMax, cost_matrix);
        }
    }

    return cost_matrix;
}


// Partition a and b (A_c and B_c in the definition) are the remining nodes to be visited, not the actual partitions
static pair<UnordSet, UnordSet> update_sets(
    UnordSet& partition_a,
    UnordSet& partition_b,
    UnordSet& current_moved_partition_a,
    UnordSet& current_moved_partition_b,
    const GainObjectImbalance& gain_object,
    const WeightedSBGraph& graph)
{
    auto node_a = UnordSet(partition_a[gain_object.i]);
    size_t partition_size_a = get_node_size(node_a, graph.get_node_weights());
    bool node_a_is_fully_used = partition_size_a == gain_object.size_i;
    if (not node_a_is_fully_used) {
        UnordSet rest_a;
        tie(node_a, rest_a) = cut_interval_by_dimension(node_a, graph.get_node_weights(), gain_object.size_i);
        cout << "cut_interval_by_dimension " << gain_object.size_i << ": " << node_a << rest_a << endl;
    }

    auto node_b = UnordSet(partition_b[gain_object.j]);
    size_t partition_size_b = get_node_size(node_b, graph.get_node_weights());
    bool node_b_is_fully_used = partition_size_b == gain_object.size_j;
    if (not node_b_is_fully_used) {
        UnordSet rest_b;
        tie(node_b, rest_b) = cut_interval_by_dimension(node_b, graph.get_node_weights(), gain_object.size_j);
        cout << "cut_interval_by_dimension " << gain_object.size_j << ": " << node_b << rest_b << endl;
    }

    partition_a = difference(partition_a, node_a);
    partition_b = difference(partition_b, node_b);

    current_moved_partition_a = cup(current_moved_partition_a, node_a);
    current_moved_partition_b = cup(current_moved_partition_b, node_b);

    flatten_set(partition_a, graph);
    flatten_set(partition_b, graph);

    flatten_set(current_moved_partition_a, graph);
    flatten_set(current_moved_partition_b, graph);

    return make_pair(node_a, node_b);
}


static void update_diff(
    CostMatrixImbalance& cost_matrix,
    UnordSet& partition_a,
    UnordSet& partition_b,
    const WeightedSBGraph& graph,
    const NodeWeight& node_weight,
    const GainObjectImbalance& gain_object,
    int LMin,
    int LMax)
{
    // this must be improved
    cost_matrix.clear();
    cost_matrix = generate_gain_matrix(graph, node_weight, partition_a, partition_b, LMin, LMax);
#if PARTITION_IMBALANCE_DEBUG
    cout << partition_a << ", " << partition_b << ", " << cost_matrix << endl;
#endif
}


int kl_sbg_imbalance(
    const WeightedSBGraph& graph,
    UnordSet& partition_a,
    UnordSet& partition_b,
    int LMin,
    int LMax)
{
#if PARTITION_IMBALANCE_DEBUG
    cout << "Algorithm starts with " << partition_a << ", " << partition_b << endl;
#endif
    auto a_c = partition_a;
    auto b_c = partition_b;
    int max_par_sum = 0;
    auto max_par_sum_set = make_pair(UnordSet(), UnordSet());
    int par_sum = 0;
    UnordSet a_v = UnordSet();
    UnordSet b_v = UnordSet();
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
        UnordSet a_, b_;
        tie(a_, b_) = update_sets(a_c, b_c, a_v, b_v, g, graph);
        update_diff(gm, a_c, b_c, graph, node_weights, g, LMin, LMax);
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
    return 0;
}


KLBipartResult kl_sbg_bipart_imbalance(const WeightedSBGraph& graph, UnordSet& partition_a,
    UnordSet& partition_b, int LMin, int LMax)
{
    auto partition_a_copy = partition_a;
    auto partition_b_copy = partition_b;
    int gain = kl_sbg_imbalance(graph, partition_a_copy, partition_b_copy, LMin, LMax);

    while (not (partition_a_copy == partition_a) and not (partition_b_copy == partition_b)) {
        partition_a = partition_a_copy;
        partition_b = partition_b_copy;
        gain = max(kl_sbg_imbalance(graph, partition_a_copy, partition_b_copy, LMin, LMax), gain);
#if PARTITION_IMBALANCE_DEBUG
        cout << "gain: " << gain << endl;
#endif
    }

#if PARTITION_IMBALANCE_DEBUG
    cout << "Final: " << partition_a << ", " << partition_b << endl;
#endif
    return KLBipartResult{partition_a, partition_b, gain};
}


}