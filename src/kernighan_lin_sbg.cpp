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

#include <set>
#include <sbg/sbg.hpp>
#include <vector>
#include <util/logger.hpp>
#include <utility>

#include <unistd.h>

#include "kernighan_lin_sbg.hpp"


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


/// This is an ad hoc function which gets the addition of the size of each interval
static unsigned get_multidim_interval_size(UnordPWMDInter& final_edges)
{
    unsigned acc = 0;
    for (size_t i = 0; i < final_edges.size(); i++) {
        for (size_t j = 0; j < final_edges[i].size(); j++) {
            acc += final_edges[i][j].end() - final_edges[i][j].begin() + 1;
        }
    }

    return acc;
}


static size_t get_c_ab(
    const UnordSet& a, const UnordSet& b,
    const BasePWMap& map_1,
    const BasePWMap& map_2)
{
    auto f = [](auto& a, auto& b, const BasePWMap& departure_map, const BasePWMap& arrival_map) {
        UnordSet comm_edges;
        for (size_t i = 0; i < departure_map.maps().size(); i++) {
            auto map_1 = *(departure_map.maps().begin() + i);
            auto map_2 = *(arrival_map.maps().begin() + i);
            auto d = preImage(a, map_1);
            auto im = image(d, map_2);
            auto inters = intersection(b, im);
            auto edges = preImage(inters, map_2);
            comm_edges = cup(edges, comm_edges);
        }

        return comm_edges;
    };

    auto intersection1 = f(a, b, map_1, map_2);
    auto intersection2 = f(a, b, map_2, map_1);

    auto communication_edges = cup(intersection1, intersection2);

    size_t comm_size = get_multidim_interval_size(communication_edges);

    return comm_size;
}


// Consider exposing these functions to facilitate being reused
ec_ic compute_EC_IC(
    const UnordSet& partition,
    const UnordSet& nodes,
    const UnordSet& partition_2,
    const BasePWMap& departure_map,
    const BasePWMap& arrival_map)
{
    UnordSet ec, ic;
    for (size_t i = 0; i < departure_map.maps().size(); i++) {
        auto map_1 = *(departure_map.maps().begin() + i);
        auto map_2 = *(arrival_map.maps().begin() + i);
        auto d = preImage(nodes, map_1);
        auto im = image(d, map_2);
        auto ic_nodes = intersection(partition, im);
        ic_nodes = difference(ic_nodes, nodes);
        auto ec_nodes = difference(im, ic_nodes);
        ec_nodes = intersection(ec_nodes, partition_2);
        auto ic_ = preImage(ic_nodes, map_2);
        auto ec_ = preImage(ec_nodes, map_2);
        ec = cup(ec_, ec);
        ic = cup(ic_, ic);
    }

    return make_pair(ec, ic);
}


static GainObject compute_diff(
    size_t i, size_t j,
    UnordSet& partition_a,
    UnordSet& partition_b,
    const BaseSBG& graph)
{
    // Firstly, copy both partitions
    auto a = partition_a[i];
    auto b = partition_b[j];

    // Take the minimum partition size. We add 1 because it includes the last element
    size_t size_node_a = a[0].end() - a[0].begin() + 1;
    size_t size_node_b = b[0].end() - b[0].begin() + 1;
    size_t min_size = min(size_node_a, size_node_b);

    // No problem here, a is just a copy of partition_a[i], same for b
    // We substract 1 because it includes the last element
    a[0].set_end(a[0].begin() + min_size - 1);
    b[0].set_end(b[0].begin() + min_size - 1);

    // Now, compute external and internal cost for both maps
    UnordSet ec_nodes_a_1, ic_nodes_a_1;
    tie(ec_nodes_a_1, ic_nodes_a_1) = compute_EC_IC(partition_a, a, partition_b, graph.map1(), graph.map2());

    UnordSet ec_nodes_a_2, ic_nodes_a_2;
    tie(ec_nodes_a_2, ic_nodes_a_2) = compute_EC_IC(partition_a, a, partition_b, graph.map2(), graph.map1());

    // Get the union between both external and internal costs for both combination of maps
    UnordSet ec_nodes_a, ic_nodes_a;
    // ec_nodes_a = intersection(cup(ec_nodes_a_1, ec_nodes_a_2), partition_b);
    ec_nodes_a = cup(ec_nodes_a_1, ec_nodes_a_2);
    ic_nodes_a = cup(ic_nodes_a_1, ic_nodes_a_2);

    // Calculate the dimension of each node
    size_t ec_a = get_multidim_interval_size(ec_nodes_a);
    size_t ic_a = get_multidim_interval_size(ic_nodes_a);
    int d_a = ec_a - ic_a;
    cout << "nodes: " << a << ": " << ec_a << ", " << ic_a << endl;

    // Same as before for partition b
    UnordSet ec_nodes_b_1, ic_nodes_b_1;
    tie(ec_nodes_b_1, ic_nodes_b_1) = compute_EC_IC(partition_b, b, partition_a, graph.map1(), graph.map2());

    UnordSet ec_nodes_b_2, ic_nodes_b_2;
    tie(ec_nodes_b_2, ic_nodes_b_2) = compute_EC_IC(partition_b, b, partition_a, graph.map2(), graph.map1());

    UnordSet ec_nodes_b, ic_nodes_b;
    // ec_nodes_b = intersection(cup(ec_nodes_b_1, ec_nodes_b_2), partition_a);
    ec_nodes_b = cup(ec_nodes_b_1, ec_nodes_b_2);
    ic_nodes_b = cup(ic_nodes_b_1, ic_nodes_b_2);

    size_t ec_b = get_multidim_interval_size(ec_nodes_b);
    size_t ic_b = get_multidim_interval_size(ic_nodes_b);
    int d_b = ec_b - ic_b;
    cout << "nodes: " << b << ": " << ec_b << ", " << ic_b << endl;

    // Get communication between a and b
    size_t c_ab = get_c_ab(a, b, graph.map1(), graph.map2());
    cout << "c_ab: " << c_ab << endl;

    // calculate gain
    int gain = d_a + d_b - 2 * c_ab;

    // Create Gain object
    const auto gain_obj = GainObject{i, j, gain, min_size};
    return gain_obj;
}



static CostMatrix generate_gain_matrix(
    UnordSet& partition_a,
    UnordSet& partition_b,
    const BaseSBG& graph)
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
    const BaseSBG& graph,
    CostMatrix& cost_matrix)
{
    // Let's add gain objects to the cost matrix
    for (size_t i = 0; i < partition_a.size(); i++) {
        auto new_gain_object = compute_diff(i, gain_object.j, partition_a, partition_b, graph);
        cost_matrix.insert(new_gain_object);
    }
}


// Partition a and b (A_c and B_c in the definition) are the remining nodes to be visited, not the actual partitions
static pair<SetPiece, SetPiece> update_sets(
    UnordSet& partition_a,
    UnordSet& partition_b,
    UnordSet& current_moved_partition_a,
    UnordSet& current_moved_partition_b,
    const GainObject& gain_object)
{
    auto node_a = partition_a[gain_object.i];
    size_t partition_size_a = node_a.intervals()[0].end() - node_a.intervals()[0].begin() + 1;
    bool node_a_is_fully_used = partition_size_a == gain_object.size;

    if (not node_a_is_fully_used) {
        node_a[0] = Interval(node_a[0].begin(), 1, node_a[0].begin() + gain_object.size - 1);
    }

    auto node_b = partition_b[gain_object.j];
    size_t partition_size_b = node_b.intervals()[0].end() - node_b.intervals()[0].begin() + 1;
    bool node_b_is_fully_used = partition_size_b == gain_object.size;
    if (not node_b_is_fully_used) {
        node_b[0] = Interval(node_b[0].begin(), 1, node_b[0].begin() + gain_object.size - 1);
    }

    // At least one of them should be fully used?
    // cout << node_a_is_fully_used << ", " <<  node_b_is_fully_used << endl;
    assert(node_a_is_fully_used or node_b_is_fully_used);

    partition_a = difference(partition_a, node_a);
    partition_b = difference(partition_b, node_b);

    current_moved_partition_a = cup(current_moved_partition_a, node_a);
    current_moved_partition_b = cup(current_moved_partition_b, node_b);

    return make_pair(node_a, node_b);
}


static GainObject max_diff(CostMatrix& cost_matrix, UnordSet& partition_a, UnordSet& partition_b, const BaseSBG& graph)
{
    // cost_matrix is sort by gain, so the first is the maximum gain
    auto g = cost_matrix.begin();

    auto gain_object = *g;
    cout << "The best is " << *g << endl;

    // remove it, we need to update those values that
    cost_matrix.erase(g);

    return gain_object;
}


static void update_diff(
    CostMatrix& cost_matrix,
    UnordSet& partition_a,
    UnordSet& partition_b,
    const BaseSBG& graph,
    const GainObject& gain_object)
{
    if (false) {
        remove_nodes_from_cost_matrix(cost_matrix, partition_a, gain_object);
        remove_nodes_from_cost_matrix(cost_matrix, partition_b, gain_object);

        update_cost_matrix(partition_b, partition_a, gain_object, graph, cost_matrix);
        update_cost_matrix(partition_a, partition_b, gain_object, graph, cost_matrix);
    }
    cost_matrix.clear();
    cost_matrix = generate_gain_matrix(partition_a, partition_b, graph);
    cout << partition_a << ", " << partition_b << ", " << cost_matrix << endl;

}


static void update_sum(
    int& par_sum,
    int g,
    int& max_par_sum,
    pair<UnordSet, UnordSet>& max_par_sum_set,
    UnordSet& a_v,
    UnordSet& b_v)
{
    par_sum += g;
    if (par_sum > max_par_sum) {
        max_par_sum = par_sum;
        max_par_sum_set = make_pair(a_v, b_v);
    }
}


void kl_sbg(const BaseSBG& graph, UnordSet& partition_a, UnordSet& partition_b)
{
    cout << "Algorithm starts with " << partition_a << ", " << partition_b << endl;
    auto a_c = partition_a;
    auto b_c = partition_b;
    int max_par_sum = 0;
    auto max_par_sum_set = make_pair(UnordSet(), UnordSet());
    int par_sum = 0;
    UnordSet a_v = UnordSet();
    UnordSet b_v = UnordSet();

    CostMatrix gm = generate_gain_matrix(partition_a, partition_b, graph);

    cout << gm << endl;

    while ((not isEmpty(a_c)) and (not isEmpty(b_c))) {
        auto g = max_diff(gm, a_c, b_c, graph);
        UnordSet a_, b_;
        tie(a_, b_) = update_sets(a_c, b_c, a_v, b_v, g);
        update_diff(gm, a_c, b_c, graph, g);
        update_sum(par_sum, g.gain, max_par_sum, max_par_sum_set, a_v, b_v);
        // sleep(2);
    }

    if (max_par_sum > 0) {
        partition_a = cup(difference(partition_a, max_par_sum_set.first), max_par_sum_set.second);
        partition_b = cup(difference(partition_b, max_par_sum_set.second), max_par_sum_set.first);
    }

    gm = generate_gain_matrix(partition_a, partition_b, graph);

    cout << "so it ends with " << gm << "\n" << max_par_sum << ", " << partition_a << ", " << partition_b << endl;

    SBG_LOG << partition_a << ", " << partition_b << endl;
}


void kl_sbg_bipart(const SBG::LIB::BaseSBG& graph, SBG::LIB::UnordSet& partition_a, SBG::LIB::UnordSet& partition_b)
{
    auto partition_a_copy = partition_a;
    auto partition_b_copy = partition_b;
    kl_sbg(graph, partition_a_copy, partition_b_copy);

    while (not (partition_a_copy == partition_a) and not (partition_b_copy == partition_b)) {
        partition_a = partition_a_copy;
        partition_b = partition_b_copy;
        kl_sbg(graph, partition_a_copy, partition_b_copy);
    }

    cout << "Final: " << partition_a << ", " << partition_b << endl;
}


}
