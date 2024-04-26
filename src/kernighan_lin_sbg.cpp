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
#include <utility>

#include "kernighan_lin_sbg.hpp"


using namespace std;

using namespace SBG::LIB;

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
       << ">";

    return os;
}


/// This is an ad hoc function which gets the addition of the size of each interval
static unsigned get_multidim_interval_size(OrdPWMDInter& final_edges)
{
    unsigned acc = 0;
    for(size_t i = 0; i < final_edges.size(); i++) {
        for(size_t j = 0; j < final_edges[i].size(); j++) {
            acc += final_edges[i][j].end() - final_edges[i][j].begin() + 1;
        }
    }

    return acc;
}


static size_t get_c_ab(
    const OrdSet& a, const OrdSet& b,
    const CanonPWMap& map_1,
    const CanonPWMap& map_2)
{
    auto f = [](auto& a, auto& b, auto& departure_map, auto& arrival_map) {
        auto d = preImage(a, departure_map);
        auto i = image(d, arrival_map);
        auto inters = intersection(b, i);
        return inters;
    };

    auto intersection1 = f(a, b, map_1, map_2);
    auto intersection2 = f(a, b, map_2, map_1);

    auto comm_nodes = cup(intersection1, intersection2);

    size_t comm_size = get_multidim_interval_size(comm_nodes);

    return comm_size;
}


ec_ic compute_EC_IC(
    const OrdSet& partition,
    const OrdSet& nodes,
    const CanonPWMap& departure_map,
    const CanonPWMap& arrival_map)
{
    auto d = preImage(nodes, departure_map);
    auto i = image(d, arrival_map);
    auto ic = intersection(partition, i);
    auto ec = difference(i, ic);

    return make_pair(ec, ic);
}


vector<ec_ic> compute_diff(
    const OrdSet& partition,
    const CanonSBG& graph)
{
    vector<ec_ic> d;
    d.reserve(partition.pieces().size());

    for (const auto& p : partition.pieces()) {
        ec_ic comunication_nodes = compute_EC_IC(partition, p, graph.map1(), graph.map2());
        ec_ic comunication_nodes_2 = compute_EC_IC(partition, p, graph.map2(), graph.map1());
        comunication_nodes.first = cup(comunication_nodes.first, comunication_nodes_2.first);
        comunication_nodes.second = cup(comunication_nodes.second, comunication_nodes_2.second);
        d.push_back(comunication_nodes);
    }

    return d;
}


static GainObject get_gain_object(
    size_t i, size_t j,
    OrdSet& partition_a,
    OrdSet& partition_b,
    const CanonSBG& graph)
{
    // Firstly, copy both partitions
    auto a = partition_a[i];
    auto b = partition_b[j];

    // Take the minimum partition size. We add 1 because it includes the last element
    size_t min_size = min(a[0].end() - a[0].begin() + 1, b[0].end() - b[0].begin() + 1);

    // No problem here, a is just a copy of partition_a[i], same for b
    // We substract 1 because it includes the last element
    a[0].set_end(a[0].begin() + min_size - 1);
    b[0].set_end(b[0].begin() + min_size - 1);

    // Now, compute external and internal cost for both maps
    OrdSet ec_nodes_a_1, ic_nodes_a_1;
    tie(ec_nodes_a_1, ic_nodes_a_1) = compute_EC_IC(partition_a, a, graph.map1(), graph.map2());

    OrdSet ec_nodes_a_2, ic_nodes_a_2;
    tie(ec_nodes_a_2, ic_nodes_a_2) = compute_EC_IC(partition_a, a, graph.map2(), graph.map1());

    // Get the union between both external and internal costs for both combination of maps
    OrdSet ec_nodes_a, ic_nodes_a;
    ec_nodes_a = cup(ec_nodes_a_1, ec_nodes_a_2);
    ic_nodes_a = cup(ic_nodes_a_1, ic_nodes_a_2);

    // Calculate the dimension of each node
    size_t ec_a = get_multidim_interval_size(ec_nodes_a);
    size_t ic_a = get_multidim_interval_size(ic_nodes_a);
    int d_a = ec_a - ic_a;

    // Same as before for partition b
    OrdSet ec_nodes_b_1, ic_nodes_b_1;
    tie(ec_nodes_b_1, ec_nodes_b_1) = compute_EC_IC(partition_b, b, graph.map1(), graph.map2());

    OrdSet ec_nodes_b_2, ic_nodes_b_2;
    tie(ec_nodes_b_2, ic_nodes_b_2) = compute_EC_IC(partition_b, b, graph.map2(), graph.map1());

    OrdSet ec_nodes_b, ic_nodes_b;
    ec_nodes_b = cup(ec_nodes_b_1, ec_nodes_b_2);
    ic_nodes_b = cup(ic_nodes_b_1, ic_nodes_b_2);

    size_t ec_b = get_multidim_interval_size(ec_nodes_b);
    size_t ic_b = get_multidim_interval_size(ic_nodes_b);
    int d_b = ec_b - ic_b;

    // Get communication between a and b
    size_t c_ab = get_c_ab(a, b, graph.map1(), graph.map2());

    // calculate gain
    int gain = d_a + d_b - 2 * c_ab;

    // Create Gain object
    const auto gain_obj = GainObject{i, j, gain, min_size};
    return gain_obj;
}


CostMatrix generate_gain_matrix(
    OrdSet& partition_a,
    OrdSet& partition_b,
    const CanonSBG& graph)
{
    CostMatrix cost_matrix;

    for (size_t i = 0; i < partition_a.pieces().size(); i++) {
        for (size_t j = 0; j < partition_b.pieces().size(); j++) {
            // Get gain object for both nodes
            auto gain_obj = get_gain_object(i, j, partition_a, partition_b, graph);

            // Insert it to the cost matrix
            cost_matrix.insert(gain_obj);
        }
    }

    return cost_matrix;
}


static bool remove_partition(size_t i, OrdSet& partition)
{
    bool was_removed = false;

    auto partition_copy = partition;
    for (size_t index = 0; index < partition.pieces().size(); index++) {
        if (index == i) {
            was_removed = true;
            continue;
        }

        partition_copy.pieces().insert(partition[index]);
    }

    partition = partition_copy;

    return was_removed;
}


static bool remove_nodes_from_cost_matrix(
    CostMatrix& cost_matrix,
    OrdSet& partition,
    const GainObject& gain_object)
{
    vector<CostMatrix::iterator> to_be_removed;
    auto node_b = partition[gain_object.j];
    size_t partition_size_b = node_b.intervals()[0].end() - node_b.intervals()[0].begin() + 1;
    bool node_b_is_fully_used = partition_size_b == gain_object.size;
    // delete each combination with partition i
    for (auto it = cost_matrix.begin(); it != cost_matrix.end(); ++it) {
        if (it->j == gain_object.j) {
            to_be_removed.push_back(it);
        }
    }

    for (auto& it : to_be_removed) {
        cout << "removing " << *it << endl;
        cost_matrix.erase(it);
    }

    return node_b_is_fully_used;
}


/// @brief update partition_b with the new interval and remove the other object from partition_a
/// @param partition_a partition to remove element with index gain_object.i
/// @param partition_b partition to update element with index gain_object.j
/// @param gain_object gain_object to be removed
static void update_partition(
    OrdSet& partition_a,
    OrdSet& partition_b,
    const GainObject& gain_object)
{
    // Firstly, create an empty partition, this one will be filled with new partition a information
    auto partition_a_copy = OrdSet();

    // take the remaining piece of node b and replace it with the previous one
    auto interval_b = partition_b[gain_object.j][0];
    partition_b[gain_object.j][0] = Interval(interval_b.begin() + gain_object.size, 1, interval_b.end());

    // Remove node from the partition
    remove_partition(gain_object.i, partition_a);

    cout << "partition a is " << partition_a << endl;
}


static void update_cost_matrix(
    OrdSet& partition_a,
    OrdSet& partition_b,
    const GainObject& gain_object,
    const CanonSBG& graph,
    CostMatrix& cost_matrix)
{
    // Let's add gain objects to the cost matrix
    for (size_t i = 0; i < partition_a.size(); i++) {
        auto new_gain_object = get_gain_object(i, gain_object.j, partition_a, partition_b, graph);
        cost_matrix.insert(new_gain_object);
    }
}


GainObject max_diff(CostMatrix& cost_matrix, OrdSet& partition_a, OrdSet& partition_b, const CanonSBG& graph)
{
    cout << "max_diff\n";
    // cost_matrix is sort by gain, so the first is the maximum gain
    auto g = cost_matrix.begin();

    auto gain_object = *g;
    cout << "The best is " << *g << endl;

    // remove it, we need to update those values that
    cost_matrix.erase(g);

    bool node_a_is_fully_used = remove_nodes_from_cost_matrix(cost_matrix, partition_a, gain_object);
    bool node_b_is_fully_used = remove_nodes_from_cost_matrix(cost_matrix, partition_b, gain_object);

    // at least one of them has to be true
    assert(node_a_is_fully_used or node_b_is_fully_used);
    // Now, check if one of the nodes is not used completelly
    // This may not be necessary but we copy partitions and update it taking
    // into account the elements to be removed.
    if (not node_a_is_fully_used) {
        update_partition(partition_b, partition_a, gain_object);
        update_cost_matrix(partition_b, partition_a, gain_object, graph, cost_matrix);
    } else {
        remove_partition(gain_object.i, partition_a);
    }

    if (not node_b_is_fully_used) {
        update_partition(partition_a, partition_b, gain_object);
        update_cost_matrix(partition_a, partition_b, gain_object, graph, cost_matrix);
    } else {
        remove_partition(gain_object.j, partition_b);
    }

    return gain_object;
}

}
