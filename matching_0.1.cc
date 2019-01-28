#include <iostream>
#include <set>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <algorithm>
using namespace std;

enum class checkStatus    { checked,unchecked };
enum class matchingStatus { matched,free };

// ---------------------------------------------
double string_to_double( const std::string& s )
// ---------------------------------------------
{
  std::istringstream i(s);
  double x;
  if (!(i >> x))
    return (double) NULL;
  return x;
}

class edge;
class node;

// ---------------------------------------------
class node
// ---------------------------------------------
{
  public:
    const basic_string<char>* namep;
    vector <edge*> incidentEdges;
    matchingStatus inM;
    node(basic_string<char>* sp, matchingStatus m) : namep(sp), inM(m) {};
    const basic_string<char> name() {
      return *namep;
    }
};

// ---------------------------------------------
class edge
// ---------------------------------------------
{
  public:
    node* a;
    node* b;
    double weight;
    matchingStatus inM;
    checkStatus inU;
    edge(node* n1, node* n2, double w, matchingStatus m, checkStatus u) : a(n1), b(n2), weight(w), inM(m), inU(u) {}
};

typedef map<string,node*> t_nodelist;
typedef vector<edge*> t_edgelist;

// ---------------------------------------------
double rand01()
// ---------------------------------------------
{
  return rand() / (RAND_MAX + 1.);
}

// ---------------------------------------------
void showNodeInfo(node* nodep)
// ---------------------------------------------
{
  cout << "  " << nodep->name();
  (nodep->inM == matchingStatus::matched) ? cout << " m" : cout << "  ";
  cout << endl;
  if (nodep->incidentEdges.size() > 0) {
    cout << "incident to edges:" << endl;
    for (unsigned i=0; i < nodep->incidentEdges.size(); i++) {
      cout << "    " << nodep->incidentEdges[i]->a->name() << "/" << nodep->incidentEdges[i]->b->name();
      (nodep->incidentEdges[i]->inM == matchingStatus::matched) ? cout << " m" : cout << "  ";
      (nodep->incidentEdges[i]->inU == checkStatus::checked) ? cout << " c" : cout << "  ";
      cout << endl;
    }
  }


}

// ---------------------------------------------
void showGraph(t_nodelist* nodesp, int bound)
// ---------------------------------------------
{
  t_nodelist nodes = *nodesp;
  map<string,node*>::iterator itr;
  itr = nodes.begin();
  unsigned c = 0;
  if (bound < 0) {
    bound = nodes.size();
  } else {
    bound = min<int>(bound,nodes.size());
  }
  while (c < bound) {
    node* nd = itr->second;
    showNodeInfo(nd);
    itr++;
    c++;
  }
  cout << endl;
}

// ---------------------------------------------
void writeMatching(t_edgelist &edges)
// ---------------------------------------------
{
  ofstream out;
  out.open ("graph_matching.txt");
  for (unsigned int i = 0; i < edges.size(); i++) {
    if (edges[i]->inM == matchingStatus::matched) { // edge is in matching
      out << edges[i]->a->name() << ' ' << edges[i]->b->name() << endl;
    }
  }
  out.close();
}

// ---------------------------------------------
void progressBar (float progress)
// ---------------------------------------------
{
  // cout << progress << " %\r";
  // cout.flush();

  int barWidth = 70;
  cout << "[";
  progress = min<double>(1,progress);
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) cout << "=";
    else if (i == pos) cout << ">";
    else cout << " ";
  }
  cout << "] " << int(progress * 100.0) << " %\r" << flush;
}

// ---------------------------------------------
edge* isEdgeUnchecked(edge* ABp,unsigned int pos)
// ---------------------------------------------
{
  node* Ap = ABp->a;
  node* Bp = ABp->b;
  t_edgelist incidentEdges;
  if (pos == 0) {
    incidentEdges = Ap->incidentEdges;
  } else {
    incidentEdges = Bp->incidentEdges;
  }
  for (unsigned int i=0; i < incidentEdges.size(); i++) {
    node* iAp = incidentEdges[i]->a;
    node* iBp = incidentEdges[i]->b;
    if ((incidentEdges[i]->inU == checkStatus::unchecked) && !((Ap == iAp) && (Bp == iBp)) && !((Ap == iBp) && (Bp == iAp))) {
      return incidentEdges[i];
    }
  }
  return NULL;
}

// ---------------------------------------------
void try_match(edge* ABp, t_edgelist &edges, double &progress)
// ---------------------------------------------
{
  t_edgelist Ca, Cb;
  float progress = 1 - (float)(unchckd.size())/edges.size();
  if (rand01() > .9) {
    progressBar(progress);
  }

  edge*      aXp   = isEdgeUnchecked(ABp,0);
  edge*      bXp   = isEdgeUnchecked(ABp,1);

  while (((ABp->a->inM == matchingStatus::free) && (ABp->b->inM == matchingStatus::free)) && ((aXp != NULL) || (bXp != NULL))) {

    if ((ABp->a->inM == matchingStatus::free) and (aXp != NULL)) {
      if (aXp->inU == checkStatus::unchecked) progress += 1;
      aXp->inU = checkStatus::checked;
      Ca.push_back(aXp);
      if (aXp->weight > ABp->weight) {
        try_match(aXp, edges, progress);
      }
    }

    if ((ABp->b->inM == matchingStatus::free) and (bXp != NULL)) {
      if (bXp->inU == checkStatus::unchecked) progress += 1;
      bXp->inU = checkStatus::checked;
      Cb.push_back(bXp);
      if (bXp->weight > ABp->weight) {
        try_match(bXp, edges, progress);
      }
    }

    aXp = isEdgeUnchecked(ABp,0);
    bXp = isEdgeUnchecked(ABp,1);

  }

  if ((ABp->a->inM == matchingStatus::matched) && (ABp->b->inM == matchingStatus::matched)) {
    // no action required (marking Ca and Cb as checkStatus::checked is already done above)
  } else if ((ABp->a->inM == matchingStatus::matched) && (ABp->b->inM == matchingStatus::free)) {
    for (unsigned int i=0; i < Cb.size(); i++) { // make edges from Cb with matchingStatus::free at end != b unchecked
      if (((Cb[i]->a == ABp->b) && (Cb[i]->b->inM == matchingStatus::free)) || ((Cb[i]->b == ABp->b) && (Cb[i]->a->inM == matchingStatus::free))) {
        Cb[i]->inU = checkStatus::unchecked;
        progress -= 1;
      }
    }
  } else if ((ABp->a->inM == matchingStatus::free) && (ABp->b->inM == matchingStatus::matched)) {
    for (unsigned int i=0; i < Ca.size(); i++) { // make edges from Ca with matchingStatus::free at end != a unchecked
      if (((Ca[i]->a == ABp->a) && (Ca[i]->b->inM == matchingStatus::free)) || ((Ca[i]->b == ABp->a) && (Ca[i]->a->inM == matchingStatus::free))) {
        Ca[i]->inU = checkStatus::unchecked;
        progress -= 1;
      }
    }
  } else /* a and b have matchingStatus::free */ {
    ABp->a->inM = matchingStatus::matched;
    ABp->b->inM = matchingStatus::matched;
    ABp->inM    = matchingStatus::matched;
    progress += 1;
  }
  ABp->inU    = checkStatus::checked;
}

// ---------------------------------------------
int pickUncheckedEdgeIndex(t_edgelist &edges, int u)
// ---------------------------------------------
{
  int i = (u+1)%edges.size();
  while (i != u) {
    if (edges[i]->inU == checkStatus::unchecked) {
      return i;
    }
    i = (i+1)%edges.size();
  }
  return -1;
}


// ---------------------------------------------
int main(int argc, char *argv[])
// ---------------------------------------------
{

  // ............................................
  srand(1);
  time_t start = time(0);
  srand (time(NULL));

  // ............................................
  t_nodelist    nodes;
  t_edgelist    edges;

  // ............................................
  cout << "Reading data... ";
  if( argc != 2 ) {
    cerr << "Please specify filename" << endl;
    return 1;
  }
  ifstream      file(argv[1]);
  if( !file ) {
    cerr << "An error occurred when opening file " << argv[1] << endl;
    return 2;
  }
  string        line;
  string        node1;
  string        node2;
  string        sweight;
  double        weight;
  t_nodelist::iterator itr;
  while(getline(file, line)) {
    stringstream  linestream(line);
    linestream >> node1 >> node2 >> sweight;
    weight = string_to_double(sweight);
    if (nodes.find(node1) == nodes.end()) { // node1 not in list
      nodes[node1] = new node(NULL,matchingStatus::free);
      itr = nodes.find(node1);
      nodes[node1]->namep = &(itr->first);
    }
    if (nodes.find(node2) == nodes.end()) { // node2 not in list
      nodes[node2] = new node(NULL,matchingStatus::free);
      t_nodelist::iterator itr = nodes.find(node2);
      nodes[node2]->namep = &(itr->first);
    }
    edge *e = new edge(nodes[node1],nodes[node2],weight,matchingStatus::free,checkStatus::unchecked);
    nodes[node1]->incidentEdges.push_back(e);
    nodes[node2]->incidentEdges.push_back(e);
    edges.push_back(e);
  }
  cout << "done." << endl;
  cout << "Nodelist contains " << nodes.size() << " nodes and " << edges.size() << " edges." << endl;

  // ............................................
  int u = pickUncheckedEdgeIndex(edges,0);
  double progress = 0;
  while (u >= 0) {  // pick edge {a,b} from U as long as there are any
    try_match(edges[u],edges,progress);
    u = pickUncheckedEdgeIndex(edges,u);
  }
  int m = 0;
  for (unsigned int i=0; i<edges.size(); i++) {
    m += (int)(edges[i]->inM == matchingStatus::matched);
  }
  cout << endl << "Matching contains " << m << " edges (=pairs) and " << (int)(m*2) << " nodes." << endl;
  cout << "We matched " << fixed << setprecision(2) << 100*m*2/nodes.size() << "% of the nodes." << endl;

  // ............................................
  cout << "writing results... ";
  writeMatching(edges);
  cout << "done." << endl;

  // ............................................
  cout << "cleaning up... ";
  nodes.clear();
  edges.clear();
  edges.shrink_to_fit();
  cout << "done." << endl;

  // ............................................
  time_t end = time(0);
  double time = difftime(end, start);
  cout << "time total: " << time << " seconds" << endl;

  return 0;
}
