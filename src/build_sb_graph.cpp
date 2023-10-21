#include <cassert>
#include <fstream>
#include <iostream>
#include <rapidjson/document.h>
#include <rapidjson/pointer.h>
#include <rapidjson/istreamwrapper.h>
#include <vector>

#include "build_sb_graph.hpp"

using namespace rapidjson;
using namespace std;

struct Var {
    string id;
    int exp_a;
    int exp_b;
    vector<int> defs;
};

struct Node {
    int id;
    int interval_start;
    int interval_end;
    vector<Var> rhs;
    vector<Var> lhs;
};


vector<Var> read_var_object(const rapidjson::Value& var_array)
{
    vector<Var> vars;
    for (const auto &value : var_array.GetArray()) {
        string id = value["id"].GetString();

        auto expression = value["exp"].GetArray();
        assert(expression.Size() == 2 && "Size of expression object is not as expected");

        continue;

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

SBG::LIB::CanonSBG build_sb_graph(const char* filename)
{
  cout << "Reading " << filename << "..." << endl;

  // Parse json document
  Document document;
  ifstream ifs(filename);
  IStreamWrapper isw(ifs);
  document.ParseStream(isw);

  // Now read the document and convert it into a known type
  auto nodes_array = document["nodes"].GetArray();
  vector<Node> nodes;
  nodes.reserve(nodes_array.Size());
  for (const auto& node : nodes_array) {
    int id = node["id"].GetInt();

    auto interval = node["interval"].GetArray();
    assert(interval.Size() == 2 && "Size of interval is not as expected");

    int interval_start = interval[0].GetInt();
    int interval_end = interval[1].GetInt();

    vector<Var> rhs = read_var_object(node["rhs"]);
    vector<Var> lhs = read_var_object(node["lhs"]);

    nodes.emplace_back(Node{id, interval_start, interval_end, rhs, lhs});
  }

  // Now, let's build a graph!
  SBG::LIB::CanonSBG graph;
  SBG::LIB::OrdSet node_set;
  for (const Node& node: nodes) {
    SBG::LIB::Interval interval(node.interval_start, 1, node.interval_end);
    node_set.emplace_hint(node_set.end(), interval);
  }

  graph.set_V(node_set);

  cout << graph << endl;

  return graph;
}