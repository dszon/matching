#include <iostream>
#include <set>
#include <string>
#include <map>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <math.h>
using namespace std;

enum class matchingStatus { matched,free };

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
    edge(node* n1, node* n2, double w, matchingStatus m) : a(n1), b(n2), weight(w), inM(m) {}};

typedef map<string,node*> t_nodelist;
typedef vector<edge*>     t_edgelist;
typedef map<edge*,edge*>  t_edgemap;




// ---------------------------------------------
double string_to_double( const string& s )
// ---------------------------------------------
{
  istringstream i(s);
  double x;
  if (!(i >> x))
    return (double) NULL;
  return x;
}

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
void showStatistics(t_edgemap &M, unsigned int nnodes)
// ---------------------------------------------
{
  int m = 0;
  double weight = 0;
  t_edgemap::iterator it;
  for (it = M.begin(); it != M.end(); it++) {
    m++;
    weight += it->second->weight;
  }
  cout << endl << "Matching contains " << m << " edges and " << (int)(m*2) << " nodes (" << fixed << setprecision(2) << 100*m*2/nnodes << "% matched). " << "Total weight: " << weight << "." << flush << endl;
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
  int barWidth = 70;
  cout << "[";
  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos) cout << "=";
    else if (i == pos) cout << ">";
    else cout << " ";
  }
  cout << "] " << ceil(100.*progress) << " %\r" << flush;
}

// ---------------------------------------------
edge* isEdgeUnchecked(edge* ABp, t_edgemap &unchckd, unsigned int pos)
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
  node* iAp;
  node* iBp;
  for (unsigned int i=0; i < incidentEdges.size(); i++) {
    iAp = incidentEdges[i]->a;
    iBp = incidentEdges[i]->b;
    if ((unchckd.find(incidentEdges[i]) != unchckd.end()) && !((Ap == iAp) && (Bp == iBp)) && !((Ap == iBp) && (Bp == iAp))) {
      return incidentEdges[i];
    }
  }
  return NULL;
}

// ---------------------------------------------
unsigned int readData(string filename,t_nodelist &nodes, t_edgelist &edges, t_edgemap &unchckd) {
// ---------------------------------------------
  ifstream      file(filename);
  if( !file ) {
    cerr << "An error occurred when opening file " << filename << endl;
    return 1;
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
    if (nodes.find(node1) == nodes.end()) {
      nodes[node1] = new node(NULL,matchingStatus::free);
      itr = nodes.find(node1);
      nodes[node1]->namep = &(itr->first);
    }
    if (nodes.find(node2) == nodes.end()) {
      nodes[node2] = new node(NULL,matchingStatus::free);
      t_nodelist::iterator itr = nodes.find(node2);
      nodes[node2]->namep = &(itr->first);
    }
    edge* ep = new edge(nodes[node1],nodes[node2],weight,matchingStatus::free);
    nodes[node1]->incidentEdges.push_back(ep);
    nodes[node2]->incidentEdges.push_back(ep);
    edges.push_back(ep);
    unchckd[ep] = ep;
  }

  cout << "done."  << flush << endl;
  cout << "Nodelist contains " << nodes.size() << " nodes and " << edges.size() << " edges."  << flush << endl;
  return 0;
}


// ---------------------------------------------
void try_match(edge* ABp, t_edgelist &edges, t_edgemap &unchckd, t_edgemap &M)
// ---------------------------------------------
{
  t_edgelist Ca, Cb;

  float progress = 1 - (float)(unchckd.size())/edges.size();

  if (rand01() > .99) {
    progressBar(progress);
  }

  edge*      aXp   = isEdgeUnchecked(ABp,unchckd,0);
  edge*      bXp   = isEdgeUnchecked(ABp,unchckd,1);

  while (((ABp->a->inM == matchingStatus::free) && (ABp->b->inM == matchingStatus::free)) && ((aXp != NULL) || (bXp != NULL))) {

    if ((ABp->a->inM == matchingStatus::free) and (aXp != NULL)) {
      unchckd.erase(aXp);
      Ca.push_back(aXp);
      if (aXp->weight > ABp->weight) {
        try_match(aXp, edges, unchckd, M);
      }
    }

    if ((ABp->b->inM == matchingStatus::free) and (bXp != NULL)) {
      unchckd.erase(bXp);
      Cb.push_back(bXp);
      if (bXp->weight > ABp->weight) {
        try_match(bXp, edges, unchckd, M);
      }
    }

    aXp = isEdgeUnchecked(ABp,unchckd,0);
    bXp = isEdgeUnchecked(ABp,unchckd,1);

  }

  if ((ABp->a->inM == matchingStatus::matched) && (ABp->b->inM == matchingStatus::matched)) {
    // no action required (marking Ca and Cb as checkStatus::checked is already done above)
  } else if ((ABp->a->inM == matchingStatus::matched) && (ABp->b->inM == matchingStatus::free)) {
    for (unsigned int i=0; i < Cb.size(); i++) { // make edges from Cb with matchingStatus::free at end != b unchecked
      if (((Cb[i]->a == ABp->b) && (Cb[i]->b->inM == matchingStatus::free)) || ((Cb[i]->b == ABp->b) && (Cb[i]->a->inM == matchingStatus::free))) {
        // Cb[i] is unchecked again
        unchckd[Cb[i]] = Cb[i];
      }
    }
  } else if ((ABp->a->inM == matchingStatus::free) && (ABp->b->inM == matchingStatus::matched)) {
    for (unsigned int i=0; i < Ca.size(); i++) { // make edges from Ca with matchingStatus::free at end != a unchecked
      if (((Ca[i]->a == ABp->a) && (Ca[i]->b->inM == matchingStatus::free)) || ((Ca[i]->b == ABp->a) && (Ca[i]->a->inM == matchingStatus::free))) {
        // Ca[i] is unchecked again
        unchckd[Ca[i]] = Ca[i];
      }
    }
  } else { // a and b have both matchingStatus::free
    ABp->a->inM = matchingStatus::matched;
    ABp->b->inM = matchingStatus::matched;
    ABp->inM    = matchingStatus::matched;
    M[ABp] = ABp;
  }
  unchckd.erase(ABp);

}

// ---------------------------------------------
void good_beta_augmentation(edge* e, t_edgelist &edges) {
// ---------------------------------------------
  node left  = e->a;
  node right = e->b;
}

// ---------------------------------------------
int main(int argc, char *argv[])
// ---------------------------------------------
{

  // doppelkanten ausschließen

  // ............................................
  srand(1);
  time_t start = time(0);
  srand (time(NULL));

  // ............................................
  t_nodelist    nodes;
  t_edgelist    edges;
  t_edgemap     unchckd, M;

  if( argc != 2 ) {
    cerr << "Please specify filename" << endl;
    return 1;
  }
  string filename = argv[1];

  // ............................................
  cout << "Reading data... " << flush;
  if (readData(filename,nodes,edges,unchckd)) {
    return 2;
  }

  // ............................................
  // Compute a 0.5-optimal solution with linear approximation scheme (Preis)
  cout << "Computing 0.5-optimal solution... " << endl;
  edge* ABp;
  while (unchckd.size() > 0) {  // pick edge {a,b} from unchcked as long as there are any
    ABp = unchckd.begin()->second; // choose first edge in list of unchecked edges to test
    try_match(ABp,edges,unchckd,M);
  }

  // ............................................
  showStatistics(M,nodes.size());

  // ............................................
  // Now we enhance the Preis-solution to get a 2/3-optimal solution (Davis and Hourgady)
  // Set M from the paper is called M here, while Mprime is realized by the inM-Variables in nodes and edges.
  cout << "Augmenting matching..." << endl;
  double beta = .8;
  t_edgemap::iterator M_it;
  for (M_it = M.begin(); M_it != M.end(); M_it++) { // visiting each edge in M exactly once
    edge* e = M_it->second;
    e->a->incidentEdges + e->b->incidentEdges;
  }


  // ............................................
  cout << "writing results... " << flush;
  writeMatching(edges);
  cout << "done."  << flush << endl;

  // ............................................
  cout << "cleaning up... " << flush;
  nodes.clear();
  edges.clear();
  edges.shrink_to_fit();
  cout << "done."  << flush << endl;

  // ............................................
  time_t end = time(0);
  double time = difftime(end, start);
  cout << "time total: " << time << " seconds" << endl;

  return 0;
}
