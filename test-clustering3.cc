#include "dna_seq.h"
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>
#include <algorithm>

using namespace std;

const int kInner = 10;
const int kOuter = 200;

char rbase() {
  int r = rand()%4;
  char b[] = "ACGT";
  return b[r];
}

string Screw(string x) {
  string ret = "";
  int cp = 0;
  for (int i = 0; i < x.size(); i++) {
    int r = rand()%100;
    if (r < 13) {
      ret += rbase();
      ret += x[i];
      cp++;
    } else if (r < 14) {
      ret += rbase();
      cp++;
    } else if (r < 15) {
    } else {
      ret += x[i];
      cp++;
    }
  }
  return ret;
}

vector<unsigned int> hh;
vector<double> jac_dist;

unsigned int KillMiddle(unsigned int x) {
  unsigned int masklow = 0xffff;
  unsigned int maskhigh = 0xffffffff << 18;
  return (x & masklow) | ((x & maskhigh) >> 2);
}

int gkm = 0;

void PrepareKmers(DNASeq& a, int in, vector<unsigned int>& ka) {
  for (int i = 0; i < a.size() - in + 1; i++) {
    ka.push_back(a.ExtractKmer(i, in));
  }
  for (int i = 0; i < a.size() - in; i++) {
    unsigned int h = KillMiddle(a.ExtractKmer(i, in+1));
    ka.push_back(h);
  }
  sort(ka.begin(), ka.end());
  auto ita = unique(ka.begin(), ka.end());
  ka.resize(distance(ka.begin(), ita));
}

bool Check(const vector<unsigned int>& ka, const vector<unsigned int>& kb, int in) {
  gkm++;
  static vector<unsigned int> inter;
  inter.clear();
  inter.resize(ka.size() + kb.size());
  auto it = set_intersection(ka.begin(), ka.end(), kb.begin(), kb.end(), inter.begin());
  int inter_size = distance(inter.begin(), it);
  double jaccard = 1.* inter_size / (ka.size() + kb.size() - inter_size);

  return jaccard >= 0.002;
}

void Test(int in, int out, string s) {
  unordered_map<unsigned int, vector<int>> index;

  vector<DNASeq> parts;
  vector<int> origin;
  vector<int> next;
  for (int c = 0; c < s.size() / 10000 * 20; c++) {
    string sc = Screw(s.substr(rand()%(s.size()-10000), 10000));
    for (int i = 0; i + out < sc.size(); i+=100) {
      parts.push_back(DNASeq(sc.substr(i, out)));
      origin.push_back(c);
      if (i > 0) {
        next.push_back(parts.size()-2);
      } else {
        next.push_back(-1);
      } 
    }
  }
  printf("parts size %d\n", parts.size());

  vector<int> order, cluster_id, rev(parts.size());
  vector<vector<unsigned int>> cluster_kmers(parts.size());
  for (int i = 0; i < parts.size(); i++) {
    order.push_back(i); cluster_id.push_back(i);
  }

  random_shuffle(order.begin(), order.end());
  vector<unsigned int> cur_kmers;

  for (int i = 0; i < order.size(); i++) {
    if (i % 10 == 0) {
      printf("i %d\r", i);
      fflush(stdout);
    }
    int el = order[i];
    rev[el] = i;
    DNASeq& ss = parts[el];
    int min_ord = i;
    cur_kmers.clear();
    PrepareKmers(ss, in, cur_kmers);
    for (int j = 0; j < ss.size() - in + 1; j++) {
      unsigned int km = ss.ExtractKmer(j, in);

      vector<int>& cands = index[km];
      for (int k = 0; k < cands.size(); k++) {
        if (cluster_id[cands[k]] != cands[k]) {
          continue;
        }
        int el_cand = order[cands[k]];
        if (origin[el] == origin[el_cand]) {
          min_ord = min(min_ord, cands[k]);
          break;
        } else if (cands[k] < min_ord) {
          if (Check(cur_kmers, cluster_kmers[cands[k]], in)) {
            min_ord = cands[k];
            break;
          }
        }
      }
    }
    for (int j = 0; j < ss.size() - in; j++) {
      unsigned int km = KillMiddle(ss.ExtractKmer(j, in + 1));

      vector<int>& cands = index[km];
      for (int k = 0; k < cands.size(); k++) {
        if (cluster_id[cands[k]] != cands[k]) {
          continue;
        }
        int el_cand = order[cands[k]];
        if (origin[el] == origin[el_cand]) {
          min_ord = min(min_ord, cands[k]);
          break;
        } else if (cands[k] < min_ord) {
          if (Check(cur_kmers, cluster_kmers[cands[k]], in)) {
            min_ord = cands[k];
            break;
          }
        }
      }
    }
    assert(cluster_id[min_ord] == min_ord);
    cluster_id[i] = min_ord;
    if (min_ord == i) {
      cluster_kmers[i] = cur_kmers;
      for (int j = 0; j < ss.size() - in + 1; j++) {
        unsigned int km = ss.ExtractKmer(j, in);
        index[km].push_back(i);
      }
      for (int j = 0; j < ss.size() - in; j++) {
        unsigned int km = KillMiddle(ss.ExtractKmer(j, in + 1));
        index[km].push_back(i);
      }
    }
  }
  unordered_map<int, int> cluster_sizes;
  int nClusters = 0;
  for (int i = 0; i < cluster_id.size(); i++) {
    if (cluster_id[i] == i) nClusters++;
    cluster_sizes[cluster_id[i]]++;
  }
  map<int, int> cluster_size_hist;
  for (auto &e: cluster_sizes) {
    cluster_size_hist[e.second]++;
  }
  printf("\n");
  for (auto &e: cluster_size_hist) {
//    printf("%d: %d\n", e.first, e.second);
  }
  printf("gkm %d\n", gkm);
  printf("clusters %d %d\n", nClusters, cluster_sizes.size());
  map<pair<int, int>, int> edges;
  set<int> verts;
  for (int i = 0; i < cluster_id.size(); i++) {
    int el = order[i];
    int ne = next[el];
    if (ne == -1) continue;
    int one = rev[ne];
    verts.insert(cluster_id[i]);
    verts.insert(cluster_id[one]);
    if (cluster_id[i] != cluster_id[one]) {
      edges[make_pair(cluster_id[i], cluster_id[one])]++;
    }
  }
  printf("verts %d\n", verts.size());
  FILE* f= fopen("wtf.dot", "w");
  fprintf(f, "digraph {\n");
  for (auto &e: verts) {
    fprintf(f, "%d [label=\"%d (%d)\"]\n",
            e, e, cluster_sizes[e]);
  }
  for (auto &e: edges) {
    fprintf(f, "%d->%d [label=\"%d\"]\n", e.first.first, e.first.second,
            e.second);
  }
  fprintf(f, "}\n");
  fclose(f);
  f= fopen("wtf2.dot", "w");
  fprintf(f, "digraph {\n");
  for (auto &e: verts) {
    fprintf(f, "%d [label=\"%d (%d)\"]\n",
            e, e, cluster_sizes[e]);
  }
  for (auto &e: edges) {
    if (e.second > 1) {
      fprintf(f, "%d->%d [label=\"%d\"]\n", e.first.first, e.first.second,
              e.second);
    }
  }
  fprintf(f, "}\n");
  fclose(f);
  
  { // compaction
    unordered_map<int, unordered_set<int>> g, gr;
    int compacted = 0;
    for (auto &e: edges) {
      if (e.second > 1 && cluster_sizes[e.first.first] >= 20 &&
          cluster_sizes[e.first.second] >= 20) {
        g[e.first.first].insert(e.first.second);
        gr[e.first.second].insert(e.first.first);
      }
    }

    printf("g size %d\n", g.size());

    // bubble removal

    bool changed = true;
    while (changed) {
      changed = false;
      vector<int> to_remove;
      for (auto &e: g) {
        for (auto &e2: e.second) {
          if (g.count(e2) == 0) continue;
          if (g[e2].size() == 1 && gr[e2].size() == 1) {
            auto next = *g[e2].begin();
            if (e.second.count(next)) {
              to_remove.push_back(e2);
              changed = true;
            }
          }
        }
      }
      for (auto &e: to_remove) {
        int prev = *gr[e].begin();
        int next = *g[e].begin();
        g[prev].erase(e);
        gr[next].erase(e);
        g.erase(e);
        gr.erase(e);
      }

      for (auto &e: g) {
        if (e.second.size() == 1) {
          auto next = *e.second.begin();
          if (next == e.first) continue;
          if (gr[next].size() == 1) {
            changed = true;
            e.second = g[next];
            g[next].clear();
            g.erase(next);
            for (auto &e2: e.second) {
              gr[e2].erase(next);
              gr[e2].insert(e.first);
            }
          }
        }
      }
    }
    printf("g size %d\n", g.size());
    FILE *of = fopen("wtf3.dot", "w");
    fprintf(of, "digraph {\n");
    for (auto &e: verts) {
      if (g.count(e) == 0) continue;
      fprintf(of, "%d [label=\"%d (%d)\"]\n",
              e, e, cluster_sizes[e]);
    }
    for (auto &e: g) {
      for (auto &e2: e.second) {
        fprintf(of, "%d -> %d\n", e.first, e2);
      }
    }
    fprintf(of, "}\n");
  }
}

int main(int argc, char** argv) {
  for (int i = 0; i < 500; i++) {
    hh.push_back(rand());
  }
  DNASeq::InitTrans();

  ifstream f(argv[1]);
  string s2, s;
  getline(f, s2);
  getline(f, s);

  Test(14, 1000, s);
}
