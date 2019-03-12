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
  cout << "Matching contains " << m << " edges and " << (int)(m*2) << " nodes (" << fixed << setprecision(2) << 100*m*2/nnodes << "% matched). " << "Total weight: " << weight << "." << flush << endl;
}

// ---------------------------------------------
void writeMatching(string filename, t_edgelist &edges)
// ---------------------------------------------
{
  ofstream out;
  out.open(filename);
  for (unsigned int i = 0; i < edges.size(); i++) {
    if (edges[i]->inM == matchingStatus::matched) { // edge is in matching
      out << edges[i]->a->name() << ',' << edges[i]->b->name() << endl;
    }
  }
  out.close();
}

// ---------------------------------------------
void progressBar (float progress)
// ---------------------------------------------
{
  // cerr << round(10000*progress)/100 << endl;
  int barWidth = 70;
  cout << "[";
  progress = progress - (int)progress;
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
  bool          deleteLoops = false;

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
  bool          alreadyThere;

  while(getline(file, line)) {
    stringstream  lineStream(line);
    getline(lineStream,node1,',');
    getline(lineStream,node2,',');
    getline(lineStream,sweight,',');

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

    // if we don't allow double edges or loops then we delete them here:
    if (deleteLoops) {
      alreadyThere = (node1 == node2);
      if (!alreadyThere) {
        for (unsigned int i=0; i < nodes[node1]->incidentEdges.size(); i++) {
          if (((nodes[node1]->incidentEdges)[i]->a->name() == node2) || ((nodes[node1]->incidentEdges)[i]->b->name() == node2)) {
            alreadyThere = true;
          }
        }
      }
    } else {
      alreadyThere = false;
    }

    if (!alreadyThere) {
      edge* ep = new edge(nodes[node1],nodes[node2],weight,matchingStatus::free);

      nodes[node1]->incidentEdges.push_back(ep);
      nodes[node2]->incidentEdges.push_back(ep);
      edges.push_back(ep);
      unchckd[ep] = ep;
    }
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

  if (rand01() > .8) {
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
  // node left  = e->a;
  // node right = e->b;
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
  t_edgemap     unchckd, M;
  double beta = 1;
  string default_sourcefile = "graph.csv";
  string sourcefile = default_sourcefile;
  string default_targetfile = "graph_matching.csv";
  string targetfile = default_targetfile;

  for (int i = 0; i < argc; ++i) {
    if (string(argv[i]) == "--in") {
      if (i+1 < argc) {
        sourcefile = argv[++i];
      } else {
        std::cerr << "ERROR: Missing filename for input (option --in)." << std::endl;
        return 1;
      }
    } else if (string(argv[i]) == "--out") {
      if (i+1 < argc) {
        targetfile = argv[++i];
      } else {
        std::cerr << "ERROR: Missing filename for output (option --out)." << std::endl;
        return 1;
      }
    } else if (string(argv[i]) == "--beta") {
      if (i+1 < argc) {
        beta = string_to_double(argv[++i]);
      } else {
        std::cerr << "ERROR: Missing value for option --beta." << std::endl;
        return 1;
      }
    }
  }

  // ............................................
  cout << "reading data ";
  if (sourcefile == default_sourcefile) {
    cout << "from default location '" << default_sourcefile << "'... " << flush;
  } else {
    cout << "from '" << sourcefile << "'... " << flush;
  }
  if (readData(sourcefile,nodes,edges,unchckd)) {
    return 1;
  }

  // ............................................
  // Compute a 0.5-optimal solution with linear approximation scheme (Preis)
  cout << "computing 0.5-optimal solution... " << endl;
  edge* ABp;
  while (unchckd.size() > 0) {  // pick edge {a,b} from unchcked as long as there are any
    ABp = unchckd.begin()->second; // choose first edge in list of unchecked edges to test
    try_match(ABp,edges,unchckd,M);
  }
  cout << endl;

  // ............................................
  showStatistics(M,nodes.size());

  // ............................................
  // Now we enhance the Preis-solution to get a 2/3-optimal solution (Davis and Hourgady)
  // Set M from the paper is called M here, while Mprime is realized by the inM-Variables in nodes and edges.
  cout << "augmenting matching with beta=" << beta << "... " << endl;

  t_edgemap::iterator M_it;
  node* elbow;
  double real_e;
  double total_gain = 0;
  vector <node*> hand(2);
  vector <edge*> best_upperarm(2);
  vector <edge*> best_lowerarm(2);
  int step = 0;

  for (M_it = M.begin(); M_it != M.end(); M_it++) { // visiting each edge in M exactly once

    float progress = 1. + (float)step/M.size();
    progressBar(progress);

    edge* e  = M_it->second;

    real_e = (int)(e->inM == matchingStatus::matched) * e->weight;
    vector <double> gain = {0,0};
    vector <double> real = {0,0};
    vector <double> alt  = {0,0};
    vector <node*> shoulders = {e->a,e->b};

    for (unsigned int s=0; s<2; s++) {
      for (unsigned int ui = 0; ui < shoulders[s]->incidentEdges.size(); ui++) {
        edge* upperarm = shoulders[s]->incidentEdges[ui];
        if ((upperarm->inM == matchingStatus::free) && (upperarm != e)) {
          alt[s]  = upperarm->weight;
          if (upperarm->a == shoulders[s]) {
            elbow = upperarm->b;
          } else {
            elbow = upperarm->a;
          }
          for (unsigned int li = 0; li < elbow->incidentEdges.size(); li++) {
            edge* lowerarm = elbow->incidentEdges[li];
            if ((lowerarm->inM == matchingStatus::matched) && (lowerarm != upperarm)) {
              real[s] = lowerarm->weight;
              if (alt[s] - beta*real[s] > gain[s]) {
                gain[s]          = alt[s] - beta*real[s];
                best_upperarm[s] = upperarm;
                best_lowerarm[s] = lowerarm;
                if (lowerarm->a == elbow) {
                  hand[s] = lowerarm->b;
                } else {
                  hand[s] = lowerarm->a;
                }
              }
            }
          }
        }
      }
    }

    if ((gain[0] + gain[1] - beta*real_e > 0) && (hand[0] != hand[1])) {
      e->inM              = matchingStatus::free;
      best_upperarm[0]->inM  = matchingStatus::matched;
      best_upperarm[1]->inM = matchingStatus::matched;
      best_lowerarm[0]->inM  = matchingStatus::free;
      best_lowerarm[1]->inM = matchingStatus::free;
      hand[0]->inM = matchingStatus::free;
      hand[1]->inM = matchingStatus::free;
      total_gain += gain[0] + gain[1] - real_e;
    }

    step++;
  }

  cout << endl << "total gain: " << total_gain   << endl;

  // ............................................
  cout << "writing results " << flush;
  if (targetfile == default_targetfile) {
    cout << "to default location '" << default_targetfile << "'... " << flush;
  } else {
    cout << "to '" << targetfile << "'... " << flush;
  }
  writeMatching(targetfile,edges);
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
  cout << "Wall time: " << time << " seconds" << endl;

  return 0;
}
