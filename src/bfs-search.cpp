#include "bfs-matcher.h"
#include "bfs-vf2.h"
#include <memory>
#include <numeric>
#include <unordered_set>

namespace glyenzy_bfs {
using glyenzy_bfs_matching::Dictionary;
using glyenzy_bfs_matching::Profile;

int root(const Profile &g) {
  return std::find(g.in.begin(), g.in.end(), 0) - g.in.begin();
}

struct Graph : Profile {
  bool alditol;
  explicit Graph(List x) : Profile(x), alditol(false) {
    List attrs = x["attributes"];
    if (attrs.containsElementNamed("alditol"))
      alditol = as<bool>(attrs["alditol"]);
  }
};

List record(const Graph &g) {
  return List::create(
      _["n"] = g.n, _["edges"] = g.edges,
      _["attributes"] = List::create(_["anomer"] = g.anomers[root(g)],
                                     _["alditol"] = g.alditol),
      _["vertices"] = List::create(_["mono"] = g.mono, _["sub"] = g.sub),
      _["edge_attributes"] = List::create(_["linkage"] = g.links));
}

struct Node {
  std::shared_ptr<Graph> compact;
  RObject graph, structure;
  string key;
  explicit Node(List x)
      : graph(x["graph"]), structure(x["structure"]),
        key(as<string>(x["key"])) {
    SEXP r = x["record"];
    if (!Rf_isNull(r))
      compact = std::make_shared<Graph>(List(r));
  }
  Node(const Graph &g, const string &k)
      : compact(std::make_shared<Graph>(g)), graph(R_NilValue),
        structure(R_NilValue), key(k) {}
  List pack() const {
    return List::create(
        _["record"] = compact ? RObject(record(*compact)) : RObject(R_NilValue),
        _["graph"] = graph, _["structure"] = structure, _["key"] = key);
  }
};

string token_position(const string &token) {
  size_t i = 0;
  while (
      i < token.size() &&
      (std::isdigit(static_cast<unsigned char>(token[i])) || token[i] == '?'))
    ++i;
  return token.substr(0, i);
}

bool contained(const string &available, const string &required, bool lenient) {
  auto a = glyenzy_bfs_matching::split(available, ','),
       r = glyenzy_bfs_matching::split(required, ',');
  if (r.size() > a.size())
    return false;
  vector<bool> used(a.size(), false);
  std::function<bool(size_t)> assign = [&](size_t i) {
    if (i == r.size())
      return true;
    for (size_t j = 0; j < a.size(); ++j) {
      string ap = token_position(a[j]), rp = token_position(r[i]);
      bool same =
          a[j] == r[i] ||
          (lenient && !ap.empty() && !rp.empty() && (ap == "?" || rp == "?") &&
           a[j].substr(ap.size()) == r[i].substr(rp.size()));
      if (!used[j] && same) {
        used[j] = true;
        if (assign(i + 1))
          return true;
        used[j] = false;
      }
    }
    return false;
  };
  return assign(0);
}

// Requirements and pruning use substituent-subset semantics, whereas enzyme
// acceptors/rejects and whole-target matching use glymotif's strict_sub
// default.
bool subset_match(const Graph &g, const Graph &m, const Dictionary &d,
                  const string &alignment, bool lenient) {
  bool plain = std::all_of(g.sub.begin(), g.sub.end(),
                           [](const string &x) { return x.empty(); }) &&
               std::all_of(m.sub.begin(), m.sub.end(),
                           [](const string &x) { return x.empty(); });
  if (plain)
    return glyenzy_bfs_matching::match(g, m, d, alignment, false, true, lenient,
                                       R_NilValue, true)
               .size() > 0;
  Graph gb = g, mb = m;
  std::fill(gb.sub.begin(), gb.sub.end(), "");
  std::fill(mb.sub.begin(), mb.sub.end(), "");
  List matches = glyenzy_bfs_matching::match(gb, mb, d, alignment, false, true,
                                             lenient, R_NilValue, false);
  for (SEXP x : matches) {
    IntegerVector mapping(x);
    bool ok = true;
    for (int i = 0; i < m.n && ok; ++i)
      ok = contained(g.sub[mapping[i] - 1], m.sub[i], lenient);
    if (ok)
      return true;
  }
  return false;
}

string label(const Graph &g, int v) {
  string sub = g.sub[v];
  sub.erase(std::remove(sub.begin(), sub.end(), ','), sub.end());
  return g.mono[v] + sub;
}

// Structural cache keys do not depend on vertex numbering or locale. Length
// prefixes prevent user residue/substituent strings from colliding with syntax.
string sized(const string &x) { return std::to_string(x.size()) + ":" + x; }
string tree_key(const Graph &g, int v) {
  vector<string> children;
  for (int e = 0; e < g.edges.nrow(); ++e)
    if (g.edges(e, 0) == v + 1)
      children.push_back(sized(g.links[e]) + tree_key(g, g.edges(e, 1) - 1));
  std::sort(children.begin(), children.end());
  string out = sized(g.mono[v]) + sized(g.sub[v]);
  for (const auto &child : children)
    out += sized(child);
  return sized(out);
}
string signature(const Graph &g) {
  return sized(g.anomers[root(g)]) + (g.alditol ? "1" : "0") +
         tree_key(g, root(g));
}

vector<int> lexical_order(const vector<string> &values, Function &order,
                          bool byte_order) {
  vector<int> out(values.size());
  std::iota(out.begin(), out.end(), 0);
  if (values.size() < 2)
    return out;
  if (byte_order) {
    std::stable_sort(out.begin(), out.end(),
                     [&](int a, int b) { return values[a] < values[b]; });
  } else {
    IntegerVector indices = order(wrap(values));
    for (int i = 0; i < indices.size(); ++i)
      out[i] = indices[i] - 1;
  }
  return out;
}

// Canonicalize once per retained product, using R's active collation for ties.
// This traversal mirrors glyrepr's canonical vertex/edge and IUPAC order.
Node canonical(const Graph &g, Function &order, bool byte_order) {
  vector<vector<int>> kids(g.n);
  vector<int> incoming(g.n, -1), depth(g.n, 0), vertices, edges;
  vector<string> sig(g.n);
  for (int e = 0; e < g.edges.nrow(); ++e) {
    int v = g.edges(e, 1) - 1;
    kids[g.edges(e, 0) - 1].push_back(v);
    incoming[v] = e;
  }
  std::function<void(int)> cache = [&](int v) {
    vector<string> tokens;
    for (int c : kids[v]) {
      cache(c);
      depth[v] = std::max(depth[v], depth[c] + 1);
      tokens.push_back(g.links[incoming[c]] + "->" + sig[c]);
    }
    sig[v] = label(g, v);
    auto ids = lexical_order(tokens, order, byte_order);
    if (!ids.empty()) {
      sig[v] += "{";
      for (size_t i = 0; i < ids.size(); ++i)
        sig[v] += (i ? "," : "") + tokens[ids[i]];
      sig[v] += "}";
    }
  };
  cache(root(g));
  auto rank = [&](int v) {
    string pos = g.links[incoming[v]].substr(3);
    return pos == "?" || pos.find('/') != string::npos ? 1.0
                                                       : 1.0 / (pos[0] - '0');
  };
  std::function<string(int)> visit = [&](int v) {
    auto children = kids[v];
    string text;
    if (!children.empty()) {
      vector<string> signatures;
      for (int c : children)
        signatures.push_back(sig[c]);
      auto ids = lexical_order(signatures, order, byte_order);
      vector<int> ranks(g.n, 0);
      int r = 0;
      for (size_t i = 0; i < ids.size(); ++i) {
        if (i && signatures[ids[i]] != signatures[ids[i - 1]])
          ++r;
        ranks[children[ids[i]]] = r;
      }
      auto better = [&](int a, int b) {
        return rank(a) != rank(b) ? rank(a) > rank(b) : ranks[a] > ranks[b];
      };
      int backbone = children[0];
      for (int c : children)
        if (depth[c] > depth[backbone] ||
            (depth[c] == depth[backbone] && better(c, backbone)))
          backbone = c;
      text = visit(backbone) + "(" + g.links[incoming[backbone]] + ")";
      edges.push_back(incoming[backbone]);
      std::stable_sort(children.begin(), children.end(), better);
      for (int c : children)
        if (c != backbone) {
          text += "[" + visit(c) + "(" + g.links[incoming[c]] + ")]";
          edges.push_back(incoming[c]);
        }
    }
    vertices.push_back(v);
    return text + label(g, v);
  };
  string key = visit(root(g)) + (g.alditol ? "-ol" : "") + "(" +
               g.anomers[root(g)] + "-";
  Graph p = g;
  vector<int> inverse(g.n);
  p.edges = IntegerMatrix(g.edges.nrow(), 2);
  for (int i = 0; i < g.n; ++i) {
    int v = vertices[i];
    inverse[v] = i;
    p.mono[i] = g.mono[v];
    p.sub[i] = g.sub[v];
    p.anomers[i] = g.anomers[v];
    p.in[i] = g.in[v];
    p.out[i] = g.out[v];
  }
  for (size_t i = 0; i < edges.size(); ++i) {
    int e = edges[i];
    p.edges(i, 0) = inverse[g.edges(e, 0) - 1] + 1;
    p.edges(i, 1) = inverse[g.edges(e, 1) - 1] + 1;
    p.links[i] = g.links[e];
  }
  return Node(p, key);
}

Graph remove_node(const Graph &g, int v) {
  Graph p = g;
  --p.n;
  p.mono.erase(p.mono.begin() + v);
  p.sub.erase(p.sub.begin() + v);
  p.anomers.erase(p.anomers.begin() + v);
  p.in.assign(p.n, 0);
  p.out.assign(p.n, 0);
  p.edges = IntegerMatrix(std::max(0, p.n - 1), 2);
  p.links.clear();
  int k = 0;
  for (int e = 0; e < g.edges.nrow(); ++e) {
    int a = g.edges(e, 0) - 1, b = g.edges(e, 1) - 1;
    if (a == v || b == v)
      continue;
    a -= a > v;
    b -= b > v;
    p.edges(k, 0) = a + 1;
    p.edges(k, 1) = b + 1;
    ++k;
    ++p.out[a];
    ++p.in[b];
    p.links.push_back(g.links[e]);
  }
  return p;
}

string optional_string(List x, const char *name) {
  SEXP value = x[name];
  return Rf_isNull(value) ? "" : as<string>(value);
}
struct Rule {
  std::shared_ptr<Graph> acceptor;
  vector<Graph> rejects, requirements;
  vector<string> alignments;
  string alignment, action, mono, linkage, sulfate, product_sub;
  int site;
  bool native = true;
  explicit Rule(List x) {
    SEXP a = x["acceptor"];
    if (Rf_isNull(a))
      native = false;
    else
      acceptor = std::make_shared<Graph>(List(a));
    List rejects_r = x["rejects"], requires_r = x["requires"];
    for (SEXP m : rejects_r) {
      if (Rf_isNull(m))
        native = false;
      else
        rejects.emplace_back(List(m));
    }
    for (SEXP r : requires_r) {
      List q(r);
      SEXP m = q["motif"];
      if (Rf_isNull(m))
        native = false;
      else {
        requirements.emplace_back(List(m));
        alignments.push_back(as<string>(q["alignment"]));
      }
    }
    site = as<int>(x["site"]) - 1;
    alignment = as<string>(x["alignment"]);
    action = as<string>(x["action"]);
    mono = optional_string(x, "mono");
    linkage = optional_string(x, "linkage");
    sulfate = optional_string(x, "sulfate");
    product_sub = optional_string(x, "product_sub");
  }
};

bool position_free(const Graph &g, int site, const string &position) {
  if (position == "?" || position.find('/') != string::npos)
    return true;
  for (int e = 0; e < g.edges.nrow(); ++e)
    if (g.edges(e, 0) == site + 1 &&
        g.links[e].substr(g.links[e].find('-') + 1) == position)
      return false;
  return true;
}
bool act(Graph &p, const Graph &g, int site, const Rule &r, bool topological) {
  if (site < 0 || site >= g.n)
    stop("Invalid native enzyme action site.");
  if (r.action == "GT") {
    if (!position_free(g, site, r.linkage.substr(r.linkage.find('-') + 1)))
      return false;
    ++p.n;
    p.mono.push_back(r.mono);
    p.sub.push_back("");
    p.anomers.push_back(r.linkage.substr(0, r.linkage.find('-')));
    p.in.push_back(1);
    p.out.push_back(0);
    ++p.out[site];
    p.links.push_back(r.linkage);
    p.edges = IntegerMatrix(g.edges.nrow() + 1, 2);
    for (int e = 0; e < g.edges.nrow(); ++e) {
      p.edges(e, 0) = g.edges(e, 0);
      p.edges(e, 1) = g.edges(e, 1);
    }
    p.edges(g.edges.nrow(), 0) = site + 1;
    p.edges(g.edges.nrow(), 1) = p.n;
  } else if (r.action == "GH") {
    if (g.n <= 1 || g.out[site])
      return false;
    p = remove_node(g, site);
  } else {
    string position = r.sulfate.substr(0, r.sulfate.size() - 1);
    if (position != "?") {
      for (const auto &s : glyenzy_bfs_matching::split(g.sub[site], ','))
        if (token_position(s) == position)
          return false;
      if (g.anomers[site].substr(1) == position ||
          !position_free(g, site, position))
        return false;
    }
    p.sub[site] = r.product_sub;
  }
  if (topological) {
    std::fill(p.links.begin(), p.links.end(),
              "?"
              "?-?");
    std::fill(p.anomers.begin(), p.anomers.end(), "??");
    p.informative = false;
  } else {
    p.informative =
        p.anomers[root(p)] != "??" ||
        std::any_of(p.links.begin(), p.links.end(), [](const string &x) {
          return x != "?"
                      "?-?";
        });
  }
  return true;
}

class Engine {
  List config;
  Dictionary dictionary;
  Function expand, filter, prune_r, target_r, order, mode_r;
  vector<Rule> rules;
  vector<List> enzymes;
  vector<std::shared_ptr<Graph>> targets;
  vector<Graph> types;
  Graph ncore, pre;
  vector<string> target_keys, remaining;
  vector<bool> active;
  std::unordered_map<string, std::shared_ptr<Node>> cache;
  Environment visited, parent, parent_enzyme, parent_step;
  vector<Node> frontier;
  vector<List> edges;
  vector<string> found;
  bool scalar, has_filter, topological, whole, product_lenient, target_lenient,
      byte_order;
  double max_steps;
  int step, candidates = 0, canonicalized = 0, cache_hits = 0, callbacks = 0;

  bool promising(const Graph &p, bool lenient) {
    // Keep R's target ordering and conditions when any target needs graph
    // metadata/floating semantics; do not skip that work after a native hit.
    if (std::any_of(targets.begin(), targets.end(),
                    [](const std::shared_ptr<Graph> &t) { return !t; }))
      return as<bool>(
          prune_r(Node(p, "").pack(), lenient ? "lenient" : "strict"));
    for (const auto &t : targets)
      if (subset_match(*t, p, dictionary, "core", target_lenient))
        return true;
    if (!subset_match(p, ncore, dictionary, "substructure", lenient) ||
        !subset_match(p, pre, dictionary, "substructure", lenient))
      return false;
    Graph trimmed = p;
    bool changed = true;
    while (changed) {
      changed = false;
      for (int v = 0; v < trimmed.n; ++v) {
        const string &mono = trimmed.mono[v];
        if (!trimmed.out[v] && trimmed.sub[v].empty() &&
            (mono == "Glc" || mono == "Man" || mono == "Hex")) {
          trimmed = remove_node(trimmed, v);
          changed = true;
          break;
        }
      }
    }
    for (const auto &t : targets)
      if (subset_match(*t, trimmed, dictionary, "core", target_lenient))
        return true;
    return false;
  }

  string glycan_type(const Graph &g) {
    if (subset_match(g, types[0], dictionary, "core", true))
      return "N";
    if (subset_match(g, types[1], dictionary, "core", true) ||
        subset_match(g, types[2], dictionary, "core", true))
      return "O";
    string m = g.mono[root(g)];
    if (m == "Man" && (subset_match(g, types[3], dictionary, "core", true) ||
                       subset_match(g, types[4], dictionary, "core", true)))
      return "";
    if (m == "GalNAc" || m == "Man" || m == "Fuc" || m == "GlcNAc" ||
        m == "Xyl")
      return "O";
    if (m == "Glc" || m == "Gal")
      return "lipid/free";
    return "";
  }

  vector<Node> generate(const Node &node, int rule_id, bool lenient) {
    const auto &r = rules[rule_id];
    const Graph &g = *node.compact;
    bool eligible = r.requirements.empty();
    for (size_t i = 0; i < r.requirements.size() && !eligible; ++i)
      eligible = subset_match(g, r.requirements[i], dictionary, r.alignments[i],
                              lenient);
    if (!eligible)
      return {};
    List matches =
        glyenzy_bfs_matching::match(g, *r.acceptor, dictionary, r.alignment,
                                    false, true, lenient, R_NilValue, false);
    vector<vector<int>> rejected;
    for (const auto &m : r.rejects) {
      List rm =
          glyenzy_bfs_matching::match(g, m, dictionary, r.alignment, false,
                                      true, lenient, R_NilValue, false);
      for (SEXP x : rm)
        rejected.push_back(as<vector<int>>(x));
    }
    vector<Node> products;
    for (SEXP m : matches) {
      vector<int> mapping = as<vector<int>>(m);
      bool reject = false;
      for (const auto &rm : rejected) {
        if (std::all_of(mapping.begin(), mapping.end(), [&](int v) {
              return std::find(rm.begin(), rm.end(), v) != rm.end();
            })) {
          reject = true;
          break;
        }
      }
      if (reject)
        continue;
      Graph p = g;
      if (!act(p, g, mapping.at(r.site) - 1, r, topological))
        continue;
      ++candidates;
      string cache_key = (lenient ? "L" : "S") + signature(p);
      auto old = cache.find(cache_key);
      if (old == cache.end()) {
        std::shared_ptr<Node> prepared;
        if (promising(p, lenient)) {
          prepared = std::make_shared<Node>(canonical(p, order, byte_order));
          ++canonicalized;
        }
        old = cache.emplace(cache_key, prepared).first;
      } else
        ++cache_hits;
      if (old->second)
        products.push_back(*old->second);
    }
    return products;
  }

  vector<int> hits(const Node &node) {
    vector<int> out;
    for (size_t j = 0; j < remaining.size(); ++j) {
      if (!active[j])
        continue;
      if (!whole) {
        if (node.key == remaining[j])
          out.push_back(j);
      } else {
        int i =
            std::find(target_keys.begin(), target_keys.end(), remaining[j]) -
            target_keys.begin();
        bool match;
        if (node.compact && targets[i])
          match = glyenzy_bfs_matching::match(*node.compact, *targets[i],
                                              dictionary, "whole", false, true,
                                              true, R_NilValue, true)
                      .size() > 0;
        else
          match = as<bool>(target_r(node.pack(), i + 1));
        if (match)
          out.push_back(j);
      }
    }
    return out;
  }
  bool unfinished() const {
    return std::any_of(active.begin(), active.end(), [](bool x) { return x; });
  }

public:
  Engine(List cfg, List cb)
      : config(cfg), dictionary(as<DataFrame>(cfg["dictionary"])),
        expand(cb["expand"]), filter(cb["filter"]), prune_r(cb["prune"]),
        target_r(cb["match_target"]), order(cb["order"]), mode_r(cb["mode"]),
        ncore(as<List>(cfg["ncore"])), pre(as<List>(cfg["pre"])),
        visited(cfg["visited"]), parent(cfg["parent"]),
        parent_enzyme(cfg["parent_enzyme"]), parent_step(cfg["parent_step"]) {
    for (SEXP r : as<List>(cfg["rules"]))
      rules.emplace_back(List(r));
    for (SEXP e : as<List>(cfg["enzymes"]))
      enzymes.emplace_back(e);
    for (SEXP t : as<List>(cfg["targets"]))
      targets.push_back(Rf_isNull(t) ? nullptr
                                     : std::make_shared<Graph>(List(t)));
    for (SEXP t : as<List>(cfg["type_motifs"]))
      types.emplace_back(List(t));
    for (SEXP x : as<List>(cfg["frontier"]))
      frontier.emplace_back(List(x));
    target_keys = as<vector<string>>(cfg["target_keys"]);
    remaining = as<vector<string>>(cfg["remaining"]);
    active.assign(remaining.size(), true);
    scalar = as<bool>(cfg["scalar"]);
    has_filter = as<bool>(cfg["filter"]);
    topological = as<bool>(cfg["topological"]);
    whole = as<bool>(cfg["whole"]);
    product_lenient = as<bool>(cfg["product_lenient"]);
    target_lenient = as<bool>(cfg["target_lenient"]);
    byte_order = as<bool>(cfg["byte_order"]);
    step = as<int>(cfg["step"]);
    max_steps = as<double>(cfg["max_steps"]);
  }

  List run() {
    Node source(as<List>(config["source"]));
    for (int h : hits(source)) {
      found.push_back(source.key);
      active[h] = false;
    }
    while (!frontier.empty() && step < max_steps && unfinished()) {
      ++step;
      checkUserInterrupt();
      vector<Node> next;
      vector<int> level_hits;
      for (const auto &current : frontier) {
        bool lenient =
            scalar ? as<bool>(mode_r(current.pack())) : product_lenient;
        std::unordered_map<int, vector<Node>> rule_results;
        string type;
        bool typed = false;
        for (size_t e = 0; e < enzymes.size(); ++e) {
          checkUserInterrupt();
          List enzyme = enzymes[e];
          string name = as<string>(enzyme["name"]);
          bool standard = as<bool>(enzyme["standard"]),
               inert = as<bool>(enzyme["inert"]);
          vector<int> ids = as<vector<int>>(enzyme["rules"]);
          bool native = standard && bool(current.compact);
          for (int id : ids)
            native = native && rules[id - 1].native;
          vector<Node> products;
          RObject originals(R_NilValue);
          if (inert)
            continue;
          if (native) {
            SEXP supported = enzyme["types"];
            if (!Rf_isNull(supported)) {
              if (!typed) {
                type = glycan_type(*current.compact);
                typed = true;
              }
              auto ts = as<vector<string>>(supported);
              bool compatible = type.empty() || std::find(ts.begin(), ts.end(),
                                                          type) != ts.end();
              if (type == "lipid/free")
                compatible =
                    std::find(ts.begin(), ts.end(), "lipid") != ts.end() ||
                    std::find(ts.begin(), ts.end(), "free") != ts.end();
              if (!compatible)
                continue;
            }
            std::unordered_set<string> cell;
            for (int id : ids) {
              auto pos = rule_results.find(id);
              if (pos == rule_results.end())
                pos =
                    rule_results.emplace(id, generate(current, id - 1, lenient))
                        .first;
              for (const auto &product : pos->second)
                if (cell.insert(product.key).second)
                  products.push_back(product);
            }
          } else {
            ++callbacks;
            List result = expand(current.pack(), e + 1);
            originals = result["products"];
            for (SEXP x : as<List>(result["nodes"]))
              products.emplace_back(List(x));
          }
          LogicalVector keep(products.size(), true);
          if (has_filter && !products.empty()) {
            List packed(products.size());
            for (size_t j = 0; j < products.size(); ++j)
              packed[j] = products[j].pack();
            keep = as<LogicalVector>(filter(packed, originals));
          }
          for (size_t j = 0; j < products.size(); ++j) {
            if (!keep[j])
              continue;
            const Node &product = products[j];
            edges.push_back(List::create(_["from"] = current.key,
                                         _["to"] = product.key,
                                         _["enzyme"] = name, _["step"] = step));
            if (!visited.exists(product.key)) {
              visited.assign(product.key, true);
              parent.assign(product.key, current.key);
              parent_enzyme.assign(product.key, name);
              parent_step.assign(product.key, step);
              next.push_back(product);
            }
            for (int h : hits(product)) {
              found.push_back(product.key);
              level_hits.push_back(h);
            }
          }
        }
      }
      for (int h : level_hits)
        active[h] = false;
      frontier = std::move(next);
    }
    vector<string> missing;
    for (size_t i = 0; i < remaining.size(); ++i)
      if (active[i])
        missing.push_back(remaining[i]);
    List queue(frontier.size()), all_edges(edges.size());
    for (size_t i = 0; i < frontier.size(); ++i)
      queue[i] = frontier[i].pack();
    for (size_t i = 0; i < edges.size(); ++i)
      all_edges[i] = edges[i];
    return List::create(
        _["found_keys"] = found, _["all_edges"] = all_edges,
        _["missing_target_keys"] = missing, _["frontier"] = queue,
        _["step"] = step,
        _["stats"] = List::create(
            _["candidates"] = candidates, _["canonicalized"] = canonicalized,
            _["cache_hits"] = cache_hits, _["callback_cells"] = callbacks));
  }
};
} // namespace glyenzy_bfs

// [[Rcpp::export]]
Rcpp::List cpp_bfs_search(Rcpp::List config, Rcpp::List callbacks) {
  return glyenzy_bfs::Engine(config, callbacks).run();
}
