#pragma once

// Resolved virtual networks use the same graph representation, rule actions,
// motif matcher and canonicalization as the forward BFS. No results persist
// across calls; shared rules are evaluated once per substrate in this call.
namespace glyenzy_bfs {
int virtual_size(const Graph &g) {
  int size = g.n;
  for (const auto &s : g.sub)
    for (const auto &token : glyenzy_bfs_matching::split(s, ','))
      size += token.back() == 'S';
  return size;
}

List virtual_trim(List source, List targets, DataFrame dictionary,
                  Function order, bool byte_order) {
  Dictionary d(dictionary);
  Node start(source);
  int from_size = virtual_size(*start.compact);
  vector<List> edges;
  vector<string> endpoints, keys;
  vector<bool> found;
  List records;
  std::unordered_set<string> saved;
  auto save = [&](const Node &node) {
    if (saved.insert(node.key).second) {
      keys.push_back(node.key);
      records.push_back(record(*node.compact));
    }
  };
  save(start);
  for (SEXP target : targets) {
    Node end{List(target)};
    int steps = virtual_size(*end.compact) - from_size;
    if (steps < 0 || !subset_match(*end.compact, *start.compact, d, "core", false)) {
      endpoints.push_back("");
      found.push_back(false);
      continue;
    }
    if (!steps) {
      endpoints.push_back(start.key);
      found.push_back(true);
      continue;
    }
    save(end);
    vector<Node> frontier{end};
    std::unordered_set<string> visited{end.key};
    bool reached = false;
    for (int depth = 0; depth < steps && !frontier.empty(); ++depth) {
      checkUserInterrupt();
      vector<Node> next;
      for (const auto &current : frontier) {
        const Graph &g = *current.compact;
        int size = virtual_size(g);
        auto keep = [&](const Graph &p, const string &enzyme) {
          if (!subset_match(p, *start.compact, d, "core", false))
            return;
          Node node = size - 1 == from_size ? start : canonical(p, order, byte_order);
          edges.push_back(List::create(_["from"] = node.key, _["to"] = current.key,
                                       _["enzyme"] = enzyme, _["step"] = double(size - from_size)));
          save(node);
          if (node.key == start.key)
            reached = true;
          else if (visited.insert(node.key).second)
            next.push_back(node);
        };
        // Preserve the R frontier's leaf-then-sulfate traversal order.
        for (int v = 0; v < g.n; ++v) {
          if (g.in[v] == 0 || g.out[v] != 0)
            continue;
          bool sulfated = false;
          for (const auto &s : glyenzy_bfs_matching::split(g.sub[v], ','))
            sulfated |= s.back() == 'S';
          if (sulfated || size - 1 < from_size)
            continue;
          for (int e = 0; e < g.edges.nrow(); ++e)
            if (g.edges(e, 1) == v + 1) {
              const string &link = g.links[e];
              keep(remove_node(g, v), link.substr(0, 1) +
                   link.substr(link.find('-') + 1) + g.mono[v] + "T");
              break;
            }
        }
        for (int v = 0; v < g.n; ++v) {
          auto tokens = glyenzy_bfs_matching::split(g.sub[v], ',');
          for (size_t i = 0; i < tokens.size(); ++i) {
            if (tokens[i].back() != 'S')
              continue;
            Graph p = g;
            p.sub[v].clear();
            for (size_t j = 0; j < tokens.size(); ++j)
              if (i != j)
                p.sub[v] += (p.sub[v].empty() ? "" : ",") + tokens[j];
            string pos = tokens[i].substr(0, tokens[i].size() - 1);
            keep(p, (pos == "3" || pos == "6" ? pos : "?") + string("SulfoT"));
          }
        }
      }
      frontier = std::move(next);
    }
    endpoints.push_back(reached ? end.key : "");
    found.push_back(reached);
  }
  records.attr("names") = wrap(keys);
  return List::create(_["found"] = found, _["endpoint_keys"] = endpoints,
                      _["all_edges"] = wrap(edges), _["records"] = records);
}

vector<Graph> virtual_rule_products(const Graph &g, const Rule &r,
                                   const Dictionary &d) {
  bool eligible = r.requirements.empty();
  for (size_t i = 0; i < r.requirements.size() && !eligible; ++i)
    eligible = subset_match(g, r.requirements[i], d, r.alignments[i], false);
  if (!eligible)
    return {};
  List matches = glyenzy_bfs_matching::match(g, *r.acceptor, d, r.alignment,
                                            false, true, false, R_NilValue, false);
  if (!matches.size())
    return {};
  vector<vector<int>> rejected;
  for (const auto &m : r.rejects) {
    List rm = glyenzy_bfs_matching::match(g, m, d, r.alignment,
                                         false, true, false, R_NilValue, false);
    for (SEXP x : rm)
      rejected.push_back(as<vector<int>>(x));
  }
  vector<Graph> products;
  for (SEXP m : matches) {
    auto mapping = as<vector<int>>(m);
    bool reject = false;
    for (const auto &rm : rejected)
      if (std::all_of(mapping.begin(), mapping.end(), [&](int v) {
            return std::find(rm.begin(), rm.end(), v) != rm.end();
          })) {
        reject = true;
        break;
      }
    if (reject)
      continue;
    Graph p = g;
    if (act(p, g, mapping.at(r.site) - 1, r, false))
      products.push_back(p);
  }
  return products;
}

// Match .glycan_types(), including its structure-string beta-Man exception.
string virtual_glycan_type(const Graph &g, const vector<Graph> &types,
                          const Dictionary &d) {
  if (subset_match(g, types[0], d, "core", true))
    return "N";
  if (subset_match(g, types[1], d, "core", true) ||
      subset_match(g, types[2], d, "core", true))
    return "O";
  string m = g.mono[root(g)];
  if (m == "Man" && !g.alditol && g.anomers[root(g)][0] == 'b')
    return "";
  if (m == "GalNAc" || m == "Man" || m == "Fuc" || m == "GlcNAc" || m == "Xyl")
    return "O";
  if (m == "Glc" || m == "Gal")
    return "lipid/free";
  return "";
}

SEXP virtual_annotate(List records, IntegerMatrix edges, List plan,
                      List type_records, DataFrame dictionary) {
  Dictionary d(dictionary);
  vector<Graph> graphs, types;
  vector<Rule> rules;
  List rule_records = plan["native_rules"], enzymes = plan["enzymes"],
       ids = plan["enzyme_rule_ids"];
  for (SEXP x : rule_records) {
    rules.emplace_back(List(x));
    if (!rules.back().native)
      return R_NilValue;
  }
  for (SEXP x : records)
    graphs.emplace_back(List(x));
  for (SEXP x : type_records)
    types.emplace_back(List(x));
  vector<string> signatures;
  for (const auto &g : graphs)
    signatures.push_back(signature(g));
  vector<vector<int>> outgoing(graphs.size());
  for (int e = 0; e < edges.nrow(); ++e)
    outgoing[edges(e, 0) - 1].push_back(e);
  vector<vector<string>> annotations(edges.nrow());
  for (size_t s = 0; s < graphs.size(); ++s) {
    if (outgoing[s].empty())
      continue;
    checkUserInterrupt();
    string type = virtual_glycan_type(graphs[s], types, d);
    std::unordered_set<string> additions;
    bool removal = false, substitution = false;
    for (int edge : outgoing[s]) {
      const Graph &target = graphs[edges(edge, 1) - 1];
      removal |= target.n == graphs[s].n - 1;
      substitution |= target.n == graphs[s].n;
      if (target.n == graphs[s].n + 1) {
        std::unordered_map<string, int> counts;
        for (int i = 0; i < target.edges.nrow(); ++i)
          ++counts[sized(target.mono[target.edges(i, 1) - 1]) + sized(target.links[i])];
        for (int i = 0; i < graphs[s].edges.nrow(); ++i)
          --counts[sized(graphs[s].mono[graphs[s].edges(i, 1) - 1]) + sized(graphs[s].links[i])];
        for (const auto &entry : counts)
          if (entry.second > 0) additions.insert(entry.first);
      }
    }
    std::unordered_map<int, std::unordered_set<string>> products;
    for (int e = 0; e < enzymes.size(); ++e) {
      List enzyme = enzymes[e];
      SEXP supported = enzyme["glycan_type"];
      if (!Rf_isNull(supported) && !type.empty()) {
        auto ts = as<vector<string>>(supported);
        auto has = [&](const string &x) { return std::find(ts.begin(), ts.end(), x) != ts.end(); };
        if (type == "lipid/free" ? !has("lipid") && !has("free") : !has(type))
          continue;
      }
      auto rule_ids = as<vector<int>>(ids[e]);
      string name = as<string>(enzyme["name"]);
      for (int id : rule_ids) {
        const Rule &rule = rules[id - 1];
        if ((rule.action == "GT" && !additions.count(sized(rule.mono) + sized(rule.linkage))) ||
            (rule.action == "GH" && !removal) ||
            (rule.action == "ST" && !substitution))
          continue;
        auto result = products.find(id);
        if (result == products.end()) {
          std::unordered_set<string> keys;
          for (const auto &p : virtual_rule_products(graphs[s], rules[id - 1], d))
            keys.insert(signature(p));
          result = products.emplace(id, std::move(keys)).first;
        }
        for (int edge : outgoing[s])
          if (result->second.count(signatures[edges(edge, 1) - 1])) {
            auto &names = annotations[edge];
            if (std::find(names.begin(), names.end(), name) == names.end())
              names.push_back(name);
          }
      }
    }
  }
  List out(edges.nrow());
  for (int e = 0; e < edges.nrow(); ++e)
    out[e] = wrap(annotations[e]);
  return out;
}
} // namespace glyenzy_bfs
