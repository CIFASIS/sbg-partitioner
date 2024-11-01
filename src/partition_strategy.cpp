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

#include "build_sb_graph.hpp"
#include "partition_strategy.hpp"

#include <bits/stdc++.h>

#define DEBUG_PARTITION_STRATEGY_ENABLED 1


using namespace std;

using namespace SBG::LIB;

namespace sbg_partitioner {


PartitionStrategyGreedy::PartitionStrategyGreedy(unsigned number_of_partitions, const WeightedSBGraph graph)
    : PartitionStrategy(),
    _number_of_partitions(number_of_partitions),
    _current_partition(0),
    _node_weight(graph.get_node_weights())
{
    // get total of nodes by accumulating all interval values
    _total_of_nodes = get_node_size(graph.V(), NodeWeight());
    size_t actual_total_of_nodes = get_node_size(graph.V(), _node_weight);

    unsigned min_amount_by_partition = actual_total_of_nodes / number_of_partitions;
    unsigned surplus = actual_total_of_nodes % number_of_partitions;
    // unsigned acceptable_surplus = ceil(min_amount_by_partition * 0.05);

    cout << "PartitionStrategyGreedy::PartitionStrategyGreedy " << actual_total_of_nodes << ", " << min_amount_by_partition << ", " << surplus << endl;
    _expected_size_by_partition = min_amount_by_partition + surplus;

    for (unsigned i = 0; i < number_of_partitions; i++) {
        _current_size_by_partition[i] = 0;
    }
}


vector<unsigned> sort_keys_by_value(const map<unsigned, unsigned>& current_size_by_partition)
{
    vector<unsigned> keys;
    keys.reserve(current_size_by_partition.size());

    for (const auto& it : current_size_by_partition) {
        keys.push_back(it.first);
    }

    sort(keys.begin(), keys.end(), [&current_size_by_partition](const auto& a, const auto& b) { return current_size_by_partition.at(a) < current_size_by_partition.at(b); });

    return keys;
}


void PartitionStrategyGreedy::operator() (const SetPiece& node)
{
#if DEBUG_PARTITION_STRATEGY_ENABLED
    cout << "Adding node " << node << endl;
#endif
    UnordSet node_to_be_added = node;

    // let's just work with sizes
    map<unsigned, unsigned> size_by_partition;
    unsigned pending_node_elements = get_node_size(node, NodeWeight());
    unsigned node_weight = get_set_cost(node, _node_weight);

    auto keys_sort_by_value = sort_keys_by_value(_current_size_by_partition);
    for (const auto i : keys_sort_by_value) {
        if (_expected_size_by_partition > _current_size_by_partition[i]) { // nothing to do for now
            auto available = _expected_size_by_partition - _current_size_by_partition[i];
            unsigned elements_can_take = available / node_weight;
            elements_can_take = min(pending_node_elements, elements_can_take);
            pending_node_elements -= elements_can_take;
            size_by_partition[i] = elements_can_take * node_weight;
            _current_size_by_partition[i] += size_by_partition[i];
        }
    }

    keys_sort_by_value = sort_keys_by_value(_current_size_by_partition);
    if (pending_node_elements > 0) {
        // look for the partion that has least elements
        unsigned p = keys_sort_by_value.front();

        unsigned to_be_added = pending_node_elements * node_weight;

        size_by_partition[p] += to_be_added;
        _current_size_by_partition[p] += to_be_added;
    }

    UnordSet remaining_node = UnordSet(node);
    for (unsigned i = 0; i < _number_of_partitions; i++) {
        if (size_by_partition[i] == 0) {
            continue;
        }

        auto &p = _partitions[i];

        UnordSet node_to_be_added;
        tie(node_to_be_added, remaining_node) = cut_interval_by_dimension(remaining_node, _node_weight, size_by_partition[i]);

        for_each(node_to_be_added.begin(), node_to_be_added.end(), [&p](const auto& set_piece) { p.insert(set_piece); });
    }

}

map<unsigned, set<SetPiece>> PartitionStrategyGreedy::partitions() const
{
    return _partitions;
}




/* PartitionStrategyDistributive */


PartitionStrategyDistributive::PartitionStrategyDistributive(unsigned number_of_partitions, const WeightedSBGraph graph)
    : PartitionStrategy(),
    _number_of_partitions(number_of_partitions),
    _nodes(graph.V()),
    _node_weight(graph.get_node_weights())
{
    for (unsigned i = 0; i < _number_of_partitions; i++) {
        auto p = make_pair(i, 0);
        _current_size_by_partition.insert(p);
    }
}


static void add_surplus_sorting_by_value(const map<unsigned, unsigned>& current_size_by_partition, map<unsigned, unsigned>& size_by_partition, unsigned surplus)
{
    // Declare vector of pairs
    vector<pair<unsigned, unsigned> > current_size_by_partition_vector;

    // Copy key-value pair from Map
    // to vector of pairs
    for (auto& it : current_size_by_partition) {
        current_size_by_partition_vector.push_back(it);
    }

    // Sort using comparator function
    sort(current_size_by_partition_vector.begin(), current_size_by_partition_vector.end(), [](const pair<unsigned, unsigned>& a, const pair<unsigned, unsigned>& b) { return a.second < b.second; });

    // Print the sorted value
    for (const auto& it : current_size_by_partition_vector) {
        cout << "surplus " << it.first << ", " << size_by_partition[it.first] << endl;
        size_by_partition[it.first]++;
        surplus--;
        if (surplus == 0) {
            break;
        }
    }
}


void PartitionStrategyDistributive::operator() (const SBG::LIB::SetPiece& node)
{
#if DEBUG_PARTITION_STRATEGY_ENABLED
    cout << "Adding " << node << " distributively to partitions" << endl;
#endif
    auto s = get_node_size(node, NodeWeight());
    unsigned size_by_part = s / _number_of_partitions;
    unsigned surplus = s % _number_of_partitions;

    map<unsigned, unsigned> size_by_partition;

    for (const auto [i, p] : _current_size_by_partition) {
        size_by_partition[i] = size_by_part;
    }

    if (surplus > 0) {
        add_surplus_sorting_by_value(_current_size_by_partition, size_by_partition, surplus);
    }

    UnordSet node_to_be_added = node;
    for (const auto [i, n] : size_by_partition) {
#if DEBUG_PARTITION_STRATEGY_ENABLED
        cout << "For partition " << i << " size: " << n << endl;
#endif
        if (size_by_partition[i] == 0) {
            continue;
        }
        auto &p = _partitions[i];
        UnordSet temp_node;

        tie(temp_node, node_to_be_added) = cut_interval_by_dimension(node_to_be_added, NodeWeight(), size_by_partition[i]);
#if DEBUG_PARTITION_STRATEGY_ENABLED
        cout << "About to add " << temp_node << " to " << i << ", remaining: " << node_to_be_added << endl;
#endif
        for_each(temp_node.pieces().begin(), temp_node.pieces().end(), [&p](const SetPiece& s) { p.insert(s); });
        _current_size_by_partition[i] += get_node_size(temp_node, _node_weight);
    }
}


map<unsigned, set<SetPiece>> PartitionStrategyDistributive::partitions() const { return _partitions; }


ostream& operator<<(ostream& os, const PartitionStrategy& pgraph)
{
    for (const auto& [id, interval] : pgraph.partitions()) {
        os << id << ": ";
        for (const auto& i : interval) {
            os << i << " ";
        }
        os << endl;
    }

    return os;
}

}