// Object matching receives compact graph records from public R accessors.
// Compatibility, floating localization, VF2 search and aggregation run here.
// Reuse the same VF2 implementation as the graph-level entry points.
#include <Rcpp.h>
#include <algorithm>
#include <set>
#include <cmath>
#include <cctype>

Rcpp::List cpp_vf2_subgraph_mono(int, Rcpp::IntegerMatrix, int, Rcpp::IntegerMatrix, Rcpp::LogicalMatrix, Rcpp::LogicalMatrix, bool, bool);

#include <unordered_map>
#include <string>
#include <functional>
using namespace Rcpp;
using std::string;
using std::vector;

namespace fused {
vector<string> split(const string& s, char sep) {
  vector<string> out; size_t start=0, pos;
  while ((pos=s.find(sep,start))!=string::npos) {
    if(pos>start) out.push_back(s.substr(start,pos-start));
    start=pos+1;
  }
  if(start<s.size()) out.push_back(s.substr(start));
  return out;
}
bool token(const string& g,const string& m,bool lenient) {
  return m=="?" || (lenient && g=="?") || g==m;
}
bool anomer(const string& g,const string& m,bool lenient) {
  return token(g.substr(0,1),m.substr(0,1),lenient) &&
    token(g.substr(1),m.substr(1),lenient);
}
bool linkage(const string& g,const string& m,bool lenient) {
  if (!anomer(g.substr(0,2),m.substr(0,2),lenient)) return false;
  auto gs=split(g.substr(3),'/'), ms=split(m.substr(3),'/');
  if(std::find(ms.begin(),ms.end(),"?")!=ms.end() ||
     (lenient && std::find(gs.begin(),gs.end(),"?")!=gs.end())) return true;
  bool any=false;
  for(auto& x:gs) {
    bool found=std::find(ms.begin(),ms.end(),x)!=ms.end();
    if(!lenient && !found) return false;
    any|=found;
  }
  return lenient?any:true;
}
bool subs(const string& g,const string& m,bool strict,bool lenient) {
  if(!strict && m.empty()) return true;
  if(g.empty() || m.empty()) return g.empty() && m.empty();
  auto gs=split(g,','),ms=split(m,',');
  if(gs.size()!=ms.size()) return false;
  vector<bool> used(ms.size(),false);
  std::function<bool(size_t)> assign=[&](size_t i) {
    if(i==gs.size()) return true;
    for(size_t j=0;j<ms.size();++j) {
      if(!used[j] && token(gs[i].substr(0,1),ms[j].substr(0,1),lenient) &&
         gs[i].substr(1)==ms[j].substr(1)) {
        used[j]=true; if(assign(i+1)) return true; used[j]=false;
      }
    }
    return false;
  };
  return assign(0);
}
struct Dictionary {
  std::unordered_map<string,string> generic;
  std::set<string> generic_names,known;
  explicit Dictionary(DataFrame table) {
    CharacterVector c=table["concrete"], g=table["generic"];
    for(int i=0;i<c.size();++i) {
      string cs=as<string>(c[i]),gs=as<string>(g[i]);
      generic.emplace(cs,gs); generic_names.insert(gs);
      known.insert(cs);known.insert(gs);
    }
  }
  bool mono(const string& g,const string& m,bool lenient) const {
    if(g==m) return true;
    bool gg=generic_names.count(g),mg=generic_names.count(m);
    auto gi=generic.find(g),mi=generic.find(m);
    if(!gg && mg && gi!=generic.end()) return gi->second==m;
    return lenient && gg && !mg && mi!=generic.end() && mi->second==g;
  }
  bool residue(const string& g,const string& gs,const string& m,
               const string& ms,bool strict,bool lenient) const {
    if(mono(g,m,lenient) && subs(gs,ms,strict,lenient)) return true;
    if(!(ms.size() && ms[0]=='?') && ms.find(",?")==string::npos) return false;
    string base,sub;
    if(g=="Neu5Ac") {base="Neu";sub="5Ac";}
    else for(auto suffix:{string("NAc"),string("N")}) {
      if(g.size()>suffix.size() && g.compare(g.size()-suffix.size(),suffix.size(),suffix)==0) {
        string b=g.substr(0,g.size()-suffix.size());
        if(known.count(b)) {base=b;sub="?"+suffix;break;}
      }
    }
    return !base.empty() && mono(base,m,lenient) &&
      subs(sub+(gs.empty()?"":","+gs),ms,strict,lenient);
  }
};
struct Profile {
  int n; vector<string> mono,sub,links,anomers;
  vector<int> in,out; IntegerMatrix edges;
  bool informative=false;
  explicit Profile(SEXP obj, bool forest=false):edges(0,2) {
    List graph(obj); n=as<int>(graph["n"]);
    List ga=graph["attributes"],va=graph["vertices"],ea=graph["edge_attributes"];
    for(auto key:{"floating_parts","floating_substituents"})
      if(!forest && ga.containsElementNamed(key) && Rf_length(ga[key])>0)
        stop("Motifs cannot contain unresolved floating metadata.");
    mono=as<vector<string>>(va["mono"]);sub=as<vector<string>>(va["sub"]);
    if(n<1 || mono.size()!=size_t(n) || sub.size()!=size_t(n)) stop("Invalid vertices.");
    string root=as<string>(ga["anomer"]);
    in.assign(n,0);out.assign(n,0);anomers.assign(n,root);
    IntegerMatrix endpoints=graph["edges"];
    IntegerVector from=endpoints(_,0),to=endpoints(_,1);
    if(from.size()!=to.size() || (!forest && from.size()!=n-1)) stop("Expected a rooted tree.");
    links=ea.containsElementNamed("linkage")?as<vector<string>>(ea["linkage"]):vector<string>();
    if(links.size()!=size_t(from.size())) stop("Invalid linkages.");
    edges=IntegerMatrix(from.size(),2); informative=root!="??";
    for(int e=0;e<from.size();++e) {
      int f=from[e]-1,t=to[e]-1;
      if(f<0||t<0||f>=n||t>=n) stop("Invalid endpoint.");
      edges(e,0)=f+1;edges(e,1)=t+1;out[f]++;in[t]++;
      anomers[t]=links[e].substr(0,links[e].find('-'));
      informative|=links[e]!="??-?";
    }
    if((!forest && std::count(in.begin(),in.end(),0)!=1) || *std::max_element(in.begin(),in.end())>1)
      stop("Expected single rooted tree.");
  }
};
struct Structures {
  vector<Profile> graphs;vector<SEXP> sources;vector<int> restore;CharacterVector codes;
  explicit Structures(SEXP obj, bool forest=false) {
    List x(obj), pool=x["graphs"];
    codes=CharacterVector(x["codes"]);
    IntegerVector index=x["restore"];
    for (int i=0;i<pool.size();++i) {
      graphs.emplace_back(pool[i],forest);
      sources.push_back(pool[i]);
    }
    for (int id:index) restore.push_back(id==NA_INTEGER ? -1 : id-1);
  }
};
List match(const Profile& g,const Profile& m,const Dictionary& dict,
           const string& alignment,bool ignore,bool strict,bool lenient,
           SEXP degree,bool first) {
  if(g.n<m.n || (alignment=="whole" && g.n!=m.n)) return List();
  bool check=!ignore && m.informative;
  if(check && !lenient && !g.informative) return List();
  bool has_degree=!Rf_isNull(degree);LogicalVector deg;
  if(has_degree) {deg=LogicalVector(degree);if(deg.size()!=m.n) stop("Degree mask length.");}
  LogicalMatrix vc(m.n,g.n),ec(m.links.size(),g.links.size());
  for(int mi=0;mi<m.n;++mi) {
    bool any=false;
    for(int gi=0;gi<g.n;++gi) {
      bool ok=dict.residue(g.mono[gi],g.sub[gi],m.mono[mi],m.sub[mi],strict,lenient);
      if(has_degree && deg[mi] == TRUE) ok=ok && m.in[mi]==g.in[gi] && m.out[mi]==g.out[gi];
      if(!has_degree && alignment=="core" && m.in[mi]==0) ok=ok && g.in[gi]==0;
      if(!has_degree && alignment=="terminal" && m.out[mi]==0) ok=ok && g.out[gi]==0;
      if(check && m.in[mi]==0) ok=ok && anomer(g.anomers[gi],m.anomers[mi],lenient);
      // Necessary degree lower bounds are valid for any subgraph monomorphism.
      ok=ok && g.in[gi]>=m.in[mi] && g.out[gi]>=m.out[mi];
      vc(mi,gi)=ok;any|=ok;
    }
    if(!any) return List();
  }
  for(int mi=0;mi<ec.nrow();++mi)for(int gi=0;gi<ec.ncol();++gi)
    ec(mi,gi)=!check || linkage(g.links[gi],m.links[mi],lenient);
  return cpp_vf2_subgraph_mono(g.n,g.edges,m.n,m.edges,vc,ec,first,!first);
}
} // namespace fused

