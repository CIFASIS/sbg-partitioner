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


#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <rapidjson/document.h>
#include <rapidjson/pointer.h>
#include <rapidjson/istreamwrapper.h>
#include <vector>
#include <util/defs.hpp>

#include "build_sb_graph.hpp"


using namespace rapidjson;
using namespace std;

using namespace SBG::LIB;
using namespace SBG::Util;


struct Var {
    string id;
    INT exp_a;
    INT exp_b;
    vector<int> defs;
};


struct Node {
    int id;
    int interval_start;
    int interval_end;
    vector<Var> rhs;
    vector<Var> lhs;
};


/// This funcion takes a json array and returns a list of parsed variable objects (Var)
static vector<Var> read_var_object(const rapidjson::Value& var_array)
{
    vector<Var> vars;
    // vars.reserve(var_array.GetArray().Size());
    for (const auto &value : var_array.GetArray()) {
        string id = value["id"].GetString();

        auto expression = value["exp"].GetArray();
        assert(expression.Size() == 2 and "Size of expression object is not as expected");

        int exp_a = expression[0].GetInt();
        int exp_b = expression[1].GetInt();

        auto def_object = value["defs"].GetArray();
        vector<int> defs;
        defs.reserve(def_object.Size());
        for(const auto& def : def_object) {
            defs.push_back(def.GetInt());
        }

        Var var = Var{id, exp_a, exp_b, defs};
        vars.push_back(var);
    }

    return vars;
}


/// This funcion takes a json object and returns a list of parsed node objects (Node)
/// and its id as key/value.
static map<int, Node> create_node_objects_from_json(const Document& document)
{
  auto nodes_array = document["nodes"].GetArray();
  map<int, Node> nodes;
  int max_counter = 0;
  for (const auto& node : nodes_array) {
    int id = node["id"].GetInt();

    auto interval = node["interval"].GetArray();
    assert(interval.Size() == 2 and "Size of interval is not as expected");

    int interval_start = interval[0].GetInt();
    int interval_end = interval[1].GetInt();
    max_counter = max(interval_end, max_counter);

    vector<Var> rhs = read_var_object(node["rhs"]);
    vector<Var> lhs = read_var_object(node["lhs"]);

    Node node_element = Node{id, interval_start, interval_end, vector<Var>(), vector<Var>()};

    for (Var var : rhs) {
      node_element.rhs.push_back(var);
    }

    for (Var var : lhs) {
      node_element.lhs.push_back(var);
    }

    nodes.insert({node_element.id, node_element});
  }

  return nodes;
}


/// Creates a set of nodes, taking into accout the offset of each one to avoid collisions.
static OrdSet create_set_of_nodes(const map<int, Node>& nodes, map<int, int>& node_offsets, int& max_value)
{
  // We start to build out set of intervals from 0
  int current_max = 0;
  OrdSet node_set;
  for (const auto& [id, node] : nodes) {

    cout << "Defining interval for node " << id << " ";

    // Define the interval and add it to the node set taking into account the current offset
    // Create an offset for each equation node. We want that each equation has its
    // own domain.
    int interval_begin = current_max;
    int interval_end = (node.interval_end - node.interval_start) + current_max;
    Interval interval(interval_begin, 1, interval_end);

    node_set.emplace_hint(node_set.end(), interval);
    cout << interval << endl;

    // set this node offset, the difference between the interval and the original one
    node_offsets[id] =  interval_begin - node.interval_start;

    // udpate max value
    current_max = interval_end + 1;
  }

  // We save max value so edge domain will not collide with node domain
  max_value = current_max;

  return node_set;
}


static vector<LExp> read_left_vars(const Node& node)
{
  const auto& l_nodes = node.lhs;

  assert(l_nodes.size() > 0);

  vector<LExp> exps(l_nodes.size());
  for (const auto& l_node : l_nodes) {
    Var var = l_node;
    LExp exp = LExp(RAT(var.exp_a, 1), RAT(var.exp_b, 1));
    exps.push_back(exp);
  }

  return exps;
}


/// This function creates a map to connect to different variables
/// @param pre_image  Subset of the domain of the expression we want to connect
/// @param edge_domain  Domain of the map
/// @param var_exp  Original expression of the variable
static CanonMap create_map_interval(const OrdPWMDInter& pre_image, const Interval& edge_domain, const LExp& var_exp, int set_offset)
{
  // If the slope is 0, we just return the expression.
  if (var_exp.slope() == 0) {
    LExp map_exp = var_exp;
    INT offset = var_exp.offset().numerator();
    offset += set_offset;
    map_exp.set_offset(RAT(offset, 1));
    return CanonMap(edge_domain, map_exp);
  }

  // Slope remains the same
  LExp map_exp;
  map_exp.set_slope(var_exp.slope());

  // We start from the original offset and substract the domain offset
  INT offset = set_offset;
  offset = offset - minElem(edge_domain);

  // Then we take the minimum element of the pre-image and add the offset
  INT min_elem = minElem(pre_image)[0];
  offset += min_elem;

  // We set the offset
  map_exp.set_offset(RAT(offset, 1));

  // And create the map
  CanonMap map = CanonMap(edge_domain, map_exp);

  return map;
}


/// Ad hoc function to get pre image of an expression from its image.
/// @param image_interval  Image we want to get the pre image from
/// @param expression  Expression to get the pre image
/// @return pre image as an interval
Interval get_pre_image(const Interval& image_interval, const LExp& expression)
{
  LExp inv_exp = inverse(expression);
  Interval pre_image = image(image_interval, inv_exp);
  return pre_image;
}


/// @brief Add documentation
/// @param nodes 
/// @param node_offsets 
/// @param max_value 
/// @return 
static tuple<OrdSet, CanonPWMap, CanonPWMap> create_graph_edges(
        const std::map<int, Node>& nodes,
        const map<int, int>& node_offsets,
        int& max_value)
{
  OrdSet edge_set;  // Our set of edges
  CanonPWMap map_left_to_right;  // Map object of one of the sides
  CanonPWMap map_right_to_left; // Map object of one of the other side

  // We iterate the set of nodes read from the input file
  for (const auto& [id, node] : nodes) {
    cout << "\nLooking for connection in node " << id << endl;

    // Define the equation interval (without offsets)
    Interval interval(node.interval_start, 1, node.interval_end);

    // Now, iterate the right expresions to connect them to their definitions.
    // When we say left_sth we are talking about the variable defined on the left side,
    // Analogous for right_sth
    for (const Var &right_var : node.rhs) {
      cout << "Looking for connection for var " << right_var.id << endl;
      LExp right_exp = LExp(RAT(right_var.exp_a, 1), RAT(right_var.exp_b, 1));
      Interval right_node_image = image(interval, right_exp);

      // Definitions of this variable are on defs field. We iterate it so we
      // can map them.
      for (int i : right_var.defs) {
        cout << "Is it connected to " << i << "?" << endl;
        auto left_node = nodes.at(i);
        Interval left_interval(left_node.interval_start, 1, left_node.interval_end);
        vector<LExp> left_exps = read_left_vars(left_node);
        for (const auto& left_exp : left_exps) {
          Interval left_node_image = image(left_interval, left_exp);

          // we want to see if the intersection of the images is not empty
          Interval image_intersection = intersection(right_node_image, left_node_image);
          if (isEmpty(image_intersection)) {
            cout << "No, it is not" << endl;
            continue;
          }
          cout << "Yes, it is" << endl;

          // Edge domain will be from the current max value and will have the quantity as the intersection
          INT domain_offset = image_intersection.end() - image_intersection.begin();

          // domain for both edges
          Interval edge_domain(max_value, 1, max_value + domain_offset);
          edge_set.emplace_hint(edge_set.end(), edge_domain);

          // Update maximum value
          max_value = maxElem(edge_domain) + 1;

          // Create and save map interval to connect right var to left var
          // right_exp.set_offset(right_exp.offset() + ); // hack
          auto pre_image_right_to_left_exp = get_pre_image(image_intersection, right_exp);
          auto right_map_interval = create_map_interval(pre_image_right_to_left_exp, edge_domain, right_exp, node_offsets.at(id));
          map_right_to_left.emplace(right_map_interval);

          // Create and save map interval to connect left var to right var
          auto pre_image_left_to_right_exp = get_pre_image(image_intersection, left_exp);
          auto left_map_interval = create_map_interval(pre_image_left_to_right_exp, edge_domain, left_exp, node_offsets.at(left_node.id));
          map_left_to_right.emplace(left_map_interval);
        }
      }
    }
  }

  return {edge_set, map_left_to_right, map_right_to_left};
}


/// @brief  Add documentation
/// @param nodes 
/// @return 
static CanonSBG create_sb_graph(const std::map<int, Node>& nodes)
{
  int max_value = 0;  // We track the max value, so we avoid domain collision between edges and nodes
  map<int, int> node_offsets;

  // Now, we create our set of nodes.
  OrdSet node_set = create_set_of_nodes(nodes, node_offsets, max_value);

  // Now, let's build a graph!
  CanonSBG graph; // This will be our graph

  const auto [edge_set, map_left_to_right, map_right_to_left] = create_graph_edges(nodes, node_offsets, max_value);

  cout << "\n";

  // set node set
  graph.set_V(node_set);

  // set edge set
  graph.set_E(edge_set);

  // set maps
  graph.set_map1(map_left_to_right);
  graph.set_map2(map_right_to_left);

  return graph;
}


CanonSBG build_sb_graph(const string& filename)
{
  cout << "Reading " << filename << "..." << endl;

  // Parse json document
  Document document;
  ifstream ifs(filename);
  IStreamWrapper isw(ifs);
  document.ParseStream(isw);

  // Now read the document and convert it into a known type
  auto nodes = create_node_objects_from_json(document);

  // Now, let's get our graph
  CanonSBG graph = create_sb_graph(nodes);

  cout << graph << endl;

  return graph;
}