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
#include "partition_strategy.hpp"


#define DEBUG_PARTITION_STRATEGY_ENABLED 0


using namespace std;

using namespace SBG::LIB;

namespace sbg_partitioner {

static unsigned get_node_size(const SetPiece& node)
{
    if (node.intervals().empty()) {
        return 0;
    }

    return node.intervals().front().end() - node.intervals().front().begin() + 1;
}


static pair<SetPiece, SetPiece> cut_interval(const SetPiece &interval, int cut_value)
{
    int interval_begin = interval.intervals().front().begin();
    int interval_end = interval.intervals().front().end();
    Interval interval_1(interval_begin, 1, cut_value);
    Interval interval_2(cut_value + 1, 1, interval_end);

    SetPiece set_1;
    set_1.emplaceBack(interval_1);

    SetPiece set_2;
    set_2.emplaceBack(interval_2);

    return make_pair(set_1, set_2);
}


PartitionStrategyGreedy::PartitionStrategyGreedy(unsigned number_of_partitions, const CanonSBG graph)
    : PartitionStrategy(),
    _number_of_partitions(number_of_partitions),
    _current_partition(0)
{
    // get total of nodes by accumulating all interval values
    _total_of_nodes = 0;
    auto nodes = graph.V();
    for (auto& node : nodes) {
        int node_size = node.intervals().front().end() - node.intervals().front().begin() + 1;
        _total_of_nodes += node_size;
    }

    unsigned min_amount_by_partition = _total_of_nodes / number_of_partitions;
    unsigned surplus = _total_of_nodes % number_of_partitions;
    unsigned acceptable_surplus = ceil(min_amount_by_partition * 0.05);

    for (unsigned i = 0; i < number_of_partitions; i++) {
        unsigned local_surplus = surplus <= acceptable_surplus ? surplus : acceptable_surplus;
        surplus -= local_surplus;
        _size_by_partition[i] = min_amount_by_partition + local_surplus;
        _current_size_by_partition[i] = 0;
    }

#if DEBUG_PARTITION_STRATEGY_ENABLED
    cout << "expected size by partition" << endl;
    for (const auto [i, s] : _size_by_partition) {
        cout << i << ", " << s << endl;
    }
#endif
}


void PartitionStrategyGreedy::operator() (const SetPiece& node)
{
#if DEBUG_PARTITION_STRATEGY_ENABLED
    cout << "Adding node " << node << endl;
#endif
    SetPiece node_to_be_added = node;
    unsigned pending_node_size = get_node_size(node);

    for (size_t i = 0; i < _number_of_partitions; i++) {

#if DEBUG_PARTITION_STRATEGY_ENABLED
        cout << "checking node " << i << endl;
#endif
        // If this happens, the partition is full, we continue with the next one
        auto &p = _partitions[i];

#if DEBUG_PARTITION_STRATEGY_ENABLED
        cout << "sizes are " << _current_size_by_partition[i] << ", " << _size_by_partition[i] << endl;
#endif
        if (_current_size_by_partition[i] == _size_by_partition[i]) {
            continue;
        }

        size_t available = _size_by_partition[i] - _current_size_by_partition[i];
        if (available < pending_node_size) {
            SetPiece temp_node;
            tie(temp_node, node_to_be_added) = cut_interval(node_to_be_added, node_to_be_added.intervals().front().begin() + available - 1);
            p.insert(temp_node);
            pending_node_size = get_node_size(node_to_be_added);
            _current_size_by_partition[i] += available;
        } else {
            p.insert(node_to_be_added);
            _current_size_by_partition[i] += pending_node_size;
            break;
        }
    }
}

map<unsigned, set<SetPiece>> PartitionStrategyGreedy::partitions() const
{
    return _partitions;
}

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