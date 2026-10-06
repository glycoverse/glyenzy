// Experimental strict, unsubstituted GT/GH BFS. No R callbacks in the search.
// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp17)]]
#include "vendor/vf2.h"
#include "vendor/matcher.h"
#include <unordered_set>
#include <cstring>
using fused::Profile;
using fused::Dictionary;

std::string tree_key(const Profile& g, int v) {
  std::vector<std::string> children;
  for(int e=0;e<g.edges.nrow();++e) if(g.edges(e,0)==v+1)
    children.push_back(g.links[e]+":"+tree_key(g,g.edges(e,1)-1));
  std::sort(children.begin(),children.end());
  std::string s=g.mono[v]+"{"+g.sub[v]+"}";
  for(auto& c:children) s+="["+c+"]";
  return s;
}
int root(const Profile& g) {return std::find(g.in.begin(),g.in.end(),0)-g.in.begin();}
std::string key(const Profile& g) {return g.anomers[root(g)]+":"+tree_key(g,root(g));}
Profile remove_node(const Profile& g,int v) {
  Profile p=g; p.n--; p.mono.erase(p.mono.begin()+v);p.sub.erase(p.sub.begin()+v);
  p.anomers.erase(p.anomers.begin()+v);p.in.assign(p.n,0);p.out.assign(p.n,0);
  p.edges=IntegerMatrix(p.n-1,2);p.links.clear();int k=0;
  for(int e=0;e<g.edges.nrow();++e) {
    int a=g.edges(e,0)-1,b=g.edges(e,1)-1;
    if(a==v||b==v)continue;
    a-=a>v;b-=b>v;p.edges(k,0)=a+1;p.edges(k,1)=b+1;
    p.out[a]++;p.in[b]++;p.links.push_back(g.links[e]);k++;
  }
  return p;
}
// Mirror glyrepr's depth/linkage/signature postorder for stable BFS ties.
Profile canonical(const Profile& g) {
  std::vector<std::vector<int>> kids(g.n);std::vector<int> incoming(g.n,-1),depth(g.n,0),vertices,edges;
  std::vector<std::string> sig(g.n);
  for(int e=0;e<g.edges.nrow();++e){int v=g.edges(e,1)-1;kids[g.edges(e,0)-1].push_back(v);incoming[v]=e;}
  auto lexical=[](const std::string& a,const std::string& b){return std::strcoll(a.c_str(),b.c_str())<0;};
  std::function<void(int)> cache=[&](int v){
    std::vector<std::string> tokens;
    for(int c:kids[v]){cache(c);depth[v]=std::max(depth[v],depth[c]+1);tokens.push_back(g.links[incoming[c]]+"->"+sig[c]);}
    sig[v]=g.mono[v];std::stable_sort(tokens.begin(),tokens.end(),lexical);
    if(!tokens.empty()){sig[v]+="{";for(size_t i=0;i<tokens.size();++i)sig[v]+=(i?",":"")+tokens[i];sig[v]+="}";}
  };cache(root(g));
  auto rank=[&](int v){auto pos=g.links[incoming[v]].substr(3);return pos=="?"||pos.find('/')!=std::string::npos?1.0:1.0/(pos[0]-'0');};
  std::function<void(int)> visit=[&](int v){
    auto children=kids[v];
    if(!children.empty()){
      auto better=[&](int a,int b){return rank(a)!=rank(b)?rank(a)>rank(b):lexical(sig[b],sig[a]);};
      int backbone=children[0];for(int c:children)if(depth[c]>depth[backbone]||(depth[c]==depth[backbone]&&better(c,backbone)))backbone=c;
      visit(backbone);edges.push_back(incoming[backbone]);std::stable_sort(children.begin(),children.end(),better);
      for(int c:children)if(c!=backbone){visit(c);edges.push_back(incoming[c]);}
    }vertices.push_back(v);
  };visit(root(g));
  Profile p=g;std::vector<int> inverse(g.n);p.edges=IntegerMatrix(g.edges.nrow(),2);
  for(int i=0;i<g.n;++i){int v=vertices[i];inverse[v]=i;p.mono[i]=g.mono[v];p.sub[i]=g.sub[v];p.anomers[i]=g.anomers[v];p.in[i]=g.in[v];p.out[i]=g.out[v];}
  for(size_t i=0;i<edges.size();++i){int e=edges[i];p.edges(i,0)=inverse[g.edges(e,0)-1]+1;p.edges(i,1)=inverse[g.edges(e,1)-1]+1;p.links[i]=g.links[e];}
  return p;
}
bool have(const Profile& g,const Profile& m,const Dictionary& d,const std::string& a) {
  return fused::match(g,m,d,a,false,false,false,R_NilValue,true).size()>0;
}
struct Rule {
  Profile acceptor; std::vector<Profile> rejects,requires;
  std::vector<std::string> require_align;std::string alignment,mono,link;
  int site;bool gt;
  Rule(List x):acceptor(x["acceptor"]) {
    alignment=as<std::string>(x["alignment"]);site=as<int>(x["site"])-1;
    gt=as<bool>(x["gt"]);mono=as<std::string>(x["mono"]);link=as<std::string>(x["link"]);
    List r=x["rejects"],q=x["requires"];
    for(SEXP a:r)rejects.emplace_back(a);
    for(SEXP a:q){List z(a);requires.emplace_back(z["motif"]);require_align.push_back(as<std::string>(z["alignment"]));}
  }
};
List record(const Profile& g) {
  return List::create(_["n"]=g.n,_["edges"]=g.edges,
    _["attributes"]=List::create(_["anomer"]=g.anomers[root(g)],_["alditol"]=false),
    _["vertices"]=List::create(_["mono"]=g.mono,_["sub"]=g.sub),
    _["edge_attributes"]=List::create(_["linkage"]=g.links));
}
// [[Rcpp::export]]
List native_bfs(List source,List targets,List enzyme_rules,DataFrame dictionary,
                List ncore,List pre,int max_steps) {
  Dictionary d(dictionary);Profile core(ncore),prem(pre);
  std::vector<Profile> tg;std::unordered_set<std::string> remaining;
  for(SEXP t:targets){tg.emplace_back(t);remaining.insert(key(tg.back()));}
  std::vector<std::vector<Rule>> enzymes;
  for(SEXP e:enzyme_rules){enzymes.emplace_back();for(SEXP r:List(e))enzymes.back().emplace_back(List(r));}
  std::vector<Profile> nodes;nodes.emplace_back(source);
  std::unordered_map<std::string,int> visited;visited[key(nodes[0])]=0;
  std::unordered_map<std::string,bool> promising;
  std::vector<int> queue{0},ef,et,ee,es,found,parents{-1},pe{-1},ps{0};
  if(remaining.erase(key(nodes[0])))found.push_back(0);
  int candidates=0,pruned=0;
  for(int step=1;step<=max_steps && !queue.empty() && !remaining.empty();++step){
    checkUserInterrupt();std::vector<int> next;std::unordered_set<std::string> hits;
    for(int id:queue){Profile g=nodes[id];
      for(size_t ei=0;ei<enzymes.size();++ei){std::unordered_set<std::string> cell;
        for(auto& r:enzymes[ei]){
          bool eligible=r.requires.empty();
          for(size_t j=0;j<r.requires.size()&&!eligible;++j)eligible=have(g,r.requires[j],d,r.require_align[j]);
          if(!eligible)continue;
          List matches=fused::match(g,r.acceptor,d,r.alignment,false,false,false,R_NilValue,false);
          std::vector<std::vector<int>> rejects;
          for(auto& reject:r.rejects){List rm=fused::match(g,reject,d,r.alignment,false,false,false,R_NilValue,false);for(SEXP m:rm)rejects.push_back(as<std::vector<int>>(m));}
          for(SEXP m:matches){std::vector<int> mapping=as<std::vector<int>>(m);bool rejected=false;
            for(auto& reject:rejects){bool subset=true;for(int v:mapping)if(std::find(reject.begin(),reject.end(),v)==reject.end()){subset=false;break;}if(subset){rejected=true;break;}}
            if(rejected)continue;int site=mapping[r.site]-1;Profile p=g;
            if(r.gt){
              std::string pos=r.link.substr(r.link.find('-')+1);bool occupied=false;
              if(pos!="?"&&pos.find('/')==std::string::npos)for(int e=0;e<g.edges.nrow();++e)if(g.edges(e,0)==site+1&&g.links[e].substr(g.links[e].find('-')+1)==pos)occupied=true;
              if(occupied)continue;
              p.n++;p.mono.push_back(r.mono);p.sub.push_back("");p.anomers.push_back(r.link.substr(0,r.link.find('-')));
              p.in.push_back(1);p.out.push_back(0);p.out[site]++;p.links.push_back(r.link);
              p.edges=IntegerMatrix(g.edges.nrow()+1,2);
              std::copy(g.edges.begin(),g.edges.begin()+g.edges.nrow(),p.edges.begin());
              for(int e=0;e<g.edges.nrow();++e)p.edges(e,1)=g.edges(e,1);
              p.edges(g.edges.nrow(),0)=site+1;p.edges(g.edges.nrow(),1)=p.n;
            }else{if(g.n<=1||g.out[site])continue;p=remove_node(g,site);}
            candidates++;std::string k=key(p);bool keep=false;auto cached=promising.find(k);
            if(cached!=promising.end())keep=cached->second;
            else{
              for(auto& target:tg)if(have(target,p,d,"core")){keep=true;break;}
              if(!keep&&have(p,core,d,"substructure")&&have(p,prem,d,"substructure")){
                Profile trimmed=p;bool changed=true;
                while(changed){changed=false;for(int v=0;v<trimmed.n;++v)if(!trimmed.out[v]&&(trimmed.mono[v]=="Man"||trimmed.mono[v]=="Glc"||trimmed.mono[v]=="Hex")){trimmed=remove_node(trimmed,v);changed=true;break;}}
                for(auto& target:tg)if(have(target,trimmed,d,"core")){keep=true;break;}
              }
              promising[k]=keep;
            }
            if(!keep){pruned++;continue;}if(!cell.insert(k).second)continue;
            int pid;auto prev=visited.find(k);
            if(prev==visited.end()){pid=nodes.size();visited[k]=pid;nodes.push_back(canonical(p));parents.push_back(id);pe.push_back(ei);ps.push_back(step);next.push_back(pid);}else pid=prev->second;
            ef.push_back(id);et.push_back(pid);ee.push_back(ei);es.push_back(step);
            if(remaining.count(k)){found.push_back(pid);hits.insert(k);}
          }
        }
      }
    }
    for(auto& h:hits)remaining.erase(h);queue=next;
  }
  List graphs(nodes.size());for(size_t i=0;i<nodes.size();++i)graphs[i]=record(nodes[i]);
  return List::create(_["graphs"]=graphs,_["from"]=ef,_["to"]=et,_["enzyme"]=ee,_["step"]=es,
    _["found"]=found,_["parent"]=parents,_["parent_enzyme"]=pe,_["parent_step"]=ps,
    _["missing"]=remaining.size(),_["candidates"]=candidates,_["pruned"]=pruned);
}
