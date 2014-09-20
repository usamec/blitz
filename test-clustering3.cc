#include "dna_seq.h"
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <map>
#include <algorithm>
#include <functional>
#include <queue>
#include "graph_manipulator.h"

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
  unsigned int masklow = 0xff;
  unsigned int maskhigh = 0xffffffff << 10;
  return (x & masklow) | ((x & maskhigh) >> 2);
}

int gkm = 0;

void radix_sort2(unsigned *begin, unsigned *end) {
    unsigned *begin1 = new unsigned[end - begin];
    unsigned *end1 = begin1 + (end - begin);

    for (unsigned shift = 0; shift < 32; shift += 8) {
        size_t count[0x100] = {};
        for (unsigned *p = begin; p != end; p++)
            count[(*p >> shift) & 0xFF]++;
        unsigned *bucket[0x100], *q = begin1;
        for (int i = 0; i < 0x100; q += count[i++])
            bucket[i] = q;
        for (unsigned *p = begin; p != end; p++)
            *bucket[(*p >> shift) & 0xFF]++ = *p;
        std::swap(begin, begin1);
        std::swap(end, end1);
    }

    delete[] begin1;
}

void uniquify(vector<unsigned int>& x) {
  unsigned begin1[3000];
  unsigned *end1 = begin1 + (x.size());
  size_t count[0x100] = {};
  for (auto &p: x)
    count[p & 0xFF]++;
  unsigned *bucket[0x100], *q = begin1;
  unsigned *bucket2[0x100];
  for (int i = 0; i < 0x100; q += count[i++]) {
    bucket[i] = q;
    bucket2[i] = q;
  }

  int out_pos = 0;
  for (int i = 0; i < x.size(); i++) {
    bool found = false;
    for (auto *p = bucket2[x[i] & 0xFF]; p != bucket[x[i] & 0xFF]; ++p) {
      if (*p == x[i]) {
        found = true;
      }
    }
    if (found) continue;
    *bucket[x[i] & 0xFF]++ = x[i];
    swap(x[i], x[out_pos]);
    out_pos++;
  }

  x.resize(out_pos);
}

void PrepareKmers(DNASeq& a, int in, vector<unsigned int>& ka) {
  for (int i = 0; i < a.size() - in + 1; i++) {
    ka.push_back(a.ExtractKmer(i, in));
  }
  for (int i = 0; i < a.size() - in; i++) {
    unsigned int h = KillMiddle(a.ExtractKmer(i, in+1));
    ka.push_back(h);
  }
//  sort(ka.begin(), ka.end());
//  radix_sort2(ka.data(), ka.data() + ka.size()); 
//  auto ita = unique(ka.begin(), ka.end());
//  ka.resize(distance(ka.begin(), ita));
  uniquify(ka);
}

bool Check(const vector<unsigned int>& ka, const vector<unsigned int>& kb, int in,
           double threshold) {
  gkm++;
/*  static vector<unsigned int> inter;
  inter.clear();
  inter.resize(ka.size() + kb.size());
  auto it = set_intersection(ka.begin(), ka.end(), kb.begin(), kb.end(), inter.begin());
  int inter_size = distance(inter.begin(), it);*/
  int inter_size = 0;
  int ia = 0, ib = 0;
  while (ia < ka.size() && ib < kb.size()) {
    if (ka[ia] == kb[ib]) {
      inter_size++;
      ia++; ib++;
    } else if (ka[ia] < kb[ib]) ia++;
    else ib++;
  }
  double jaccard = 1.* inter_size / (ka.size() + kb.size() - inter_size);

  return jaccard >= threshold;
}

pair<int, int> HasSubseq(const vector<int>& subseq, const vector<int>& cluster,
                         const unordered_set<int>& active_clusters) {
  int first_match = cluster.size();
  int last = 0;
  for (int i = 0; i < subseq.size(); i++) {
    int pos = -1;
    for (int j = last; j < cluster.size(); j++) {
      if (cluster[j] == subseq[i]) {
        first_match = min(first_match, j);
        pos = j;
        break;
      }
    }
    if (pos == -1) return make_pair(-1, -1);
    last = pos+1;
  }
  while (last < cluster.size() && cluster[last] == subseq.back()) {
    last++;
  }
  while (last < cluster.size() && active_clusters.count(cluster[last]) == 0) {
    last++;
  }
  return make_pair(first_match, last);
}

string RunFalcon(vector<int>& good_reads, vector<pair<int, int>>& good_reads_pos,
               vector<string>& reads) {
  int longest = 0;
  assert(good_reads.size() == good_reads_pos.size());
  for (int i = 0; i < good_reads.size(); i++) {
    if (good_reads_pos[longest].second - good_reads_pos[longest].first <
        good_reads_pos[i].second - good_reads_pos[i].first) {
      longest = i;
    }
  }
  printf("starting falcon %d ", good_reads.size());

  FILE* f = fopen("tmp.in", "w");
  fprintf(f, "seed %s\n", reads[good_reads[longest]].substr(
          good_reads_pos[longest].first*100,
          (good_reads_pos[longest].second - good_reads_pos[longest].first) * 100).c_str());
  for (int i = 0; i < good_reads.size(); i++) {
    fprintf(f, "read%d %s\n", i, reads[good_reads[i]].c_str());/*substr(
          good_reads_pos[i].first*100,
          (good_reads_pos[i].second - good_reads_pos[i].first) * 100).c_str());*/
  }
  fprintf(f, "+ +\n- -\n");
  fclose(f);
  printf("falcon ready\n");
  system("falcon_sense.py --min_cov 2 <tmp.in >tmp.out");
  ifstream of("tmp.out");
  string s;
  if (!getline(of, s)) {
    printf("fail\n");
    return "";
  }
  if (!getline(of, s)) {
    printf("fail\n");
    return "";
  }
  printf("got %d\n", s.length());
//  printf("start %s\n", s.substr(0, 100).c_str());
//  printf("end %s\n", s.substr(s.size()-100, 100).c_str());
  return s;
}

class Falcon {
 public:
  Falcon() {
    f = fopen("tmp.in", "w");
    id = 0;
  }

  void AddSeed(vector<int>& good_reads, vector<pair<int, int>>& good_reads_pos,
               vector<string>& reads) {
    int longest = 0;
    assert(good_reads.size() == good_reads_pos.size());
    for (int i = 0; i < good_reads.size(); i++) {
      if (good_reads_pos[longest].second - good_reads_pos[longest].first <
          good_reads_pos[i].second - good_reads_pos[i].first) {
        longest = i;
      }
    }
    fprintf(f, "%d %s\n", id, reads[good_reads[longest]].substr(
            good_reads_pos[longest].first*100,
            (good_reads_pos[longest].second - good_reads_pos[longest].first) * 100).c_str());
    for (int i = 0; i < good_reads.size(); i++) {
      fprintf(f, "read%d_%d %s\n", id, i, reads[good_reads[i]].c_str());/*substr(
            good_reads_pos[i].first*100,
            (good_reads_pos[i].second - good_reads_pos[i].first) * 100).c_str());*/
    }
    fprintf(f, "+ +\n");
    id++;
  }

  vector<string> Run() {
    vector<string> ret(id);
    fprintf(f, "- -\n");
    fclose(f);
    printf("falcon ready\n");
    system("falcon_sense.py --n_core 4 --min_cov 2 <tmp.in >tmp.out");
    ifstream of("tmp.out");
    string s, s2;
    while (getline(of, s)) {
      getline(of, s2);
      int num = atoi(s.substr(1).c_str());
      ret[num] = s2;
    }
    return ret;
  }

  FILE *f;
  int id;
};

int Overlap(string x, string last) {
  vector<vector<int>> mat(x.size() + 1, vector<int>(last.size()+1));
  // mat[0][i] -- any start place in last -> 0
  // mat[i][0] -- i-th posibion in x -> -i
  for (int i = 1; i <= x.size(); i++) {
    mat[i][0] = -i;
    for (int j = 1; j <= last.size(); j++) {
      mat[i][j] = max(
        mat[i-1][j] - 1,
        max(mat[i][j-1] - 1,
            mat[i-1][j-1] - 1));
      if (x[i-1] == last[j-1]) {
        mat[i][j] = max(mat[i][j], mat[i-1][j-1] + 1);
      }
    }
  }
  // mat[i][vela] - last column in last
  int best_pos = 0;
  for (int j = 0; j <= x.size(); j++) {
    if (mat[j][last.size()] < j / 2) continue;
    if (mat[best_pos][last.size()] < mat[j][last.size()]) {
      best_pos = j;
    }
  }
  printf("best %d %d\n", best_pos, mat[best_pos][last.size()]);
  return best_pos;
}

pair<int,int> FastOverlap(string x, string& prev) {
  unordered_map<string, vector<int>> xindex;
  int k = 30;
  for (int i = 0; i + k <= x.size(); i++) {
    xindex[x.substr(i, k)].push_back(i);
  }

  for (int i = prev.size() - k; i >= prev.size() - k - 1000; i--) {
    string probe = prev.substr(i, k);
    if (xindex.count(probe) && xindex[probe].size() == 1) {
      return make_pair(xindex[probe][0] + k, i+k);
    }
  }
  return make_pair(0, 0);
}

void ProcessGraph(vector<vector<int>>& read_clusters, 
                  map<pair<int, int>, int>& edges,
                  unordered_map<int, int>& cluster_sizes, 
                  vector<vector<unsigned int>>& cluster_kmers,
                  vector<string>& reads) {
  GraphManipulator gm(read_clusters, edges, cluster_sizes, cluster_kmers);
  gm.OutputGraph("wtf.dot");
  gm.Process();
  gm.OutputGraph("wtf2.dot", false);

  FILE *of = fopen("of.fasta", "w");
  FILE *of2 = fopen("of2.fasta", "w");
  unordered_set<int> active_clusters;
  for (auto &e: gm.in_node_seq) {
    active_clusters.insert(e.begin(), e.end());
  }
  for (int i = 0; i < gm.in_node_seq.size(); i++) {
    vector<int>& seq = gm.in_node_seq[i];
    reverse(seq.begin(), seq.end());
    if (seq.size() < 5) continue;
    string last = "";
    int pos = 0;
    Falcon falcon;
    for (int j = 0; j + 1 < seq.size(); j++) {
//      printf("i %d j %d %d %d\n", i, j, seq[j], pos);
      vector<int> good_reads;
      vector<pair<int, int>> good_reads_pos;
      for (int k = 0; k < read_clusters.size(); k++) {
        pair<int, int> subseq_pos = HasSubseq(
            vector<int>(seq.begin()+j, seq.begin()+(j+2)), read_clusters[k],
            active_clusters);
        if (subseq_pos.first != -1) {
          good_reads.push_back(k);
          good_reads_pos.push_back(subseq_pos);
        }
      }
      if (good_reads.size() == 0) {
        printf("wat\n");
      }
      falcon.AddSeed(good_reads, good_reads_pos, reads);
/*      string x = RunFalcon(good_reads, good_reads_pos, reads);
      if (x.size() > 0) {
        fprintf(of, ">%d_%d\n", i, j);
        fprintf(of, "%s\n", x.c_str());
      }*/
/*      if (j == 0) {
        fprintf(of2, ">output_%d\n", i);
      }
      if (last.size() > 0) {
        int move = Overlap(x, last);
        fprintf(of2, "%s\n", x.substr(move).c_str());
        pos += x.size() - move;
      } else {
        fprintf(of2, "%s\n", x.c_str());
        pos += x.size();
      }*/
      //last = x;
    }
    vector<string> parts = falcon.Run();
    string output;
    for (int j = 0; j < parts.size(); j++) {
      if (parts[j].size() > 0) {
        fprintf(of, ">%d_%d\n", i, j);
        fprintf(of, "%s\n", parts[j].c_str());
        if (j == 0) {
          output += parts[j];
        } else {
          int back_move = min(output.size(), parts[j].size()+1000);
//          int move = Overlap(parts[j], output.substr(output.size() - back_move));
          pair<int,int> move2 = FastOverlap(parts[j], output);
//          printf("move %d %d\n", move, move2.first);
          if (move2.first < 50) {
            output += string(1000, 'N');
            output += parts[j];
          } else {
            output.resize(move2.second);
            output += parts[j].substr(move2.first);
          }
        }
      } else {
        output += string(1000, 'N');
      }
    }
    fprintf(of2, ">%d\n", i);
    fprintf(of2, "%s\n", output.c_str());
  }
  fclose(of);
  fclose(of2);
}

void Test(int in, int out, string s) {
  unordered_map<unsigned int, vector<int>> index;

  vector<DNASeq> parts;
  vector<int> origin;
  vector<int> origin_order;
  vector<int> next;
  vector<vector<int>> read_clusters;
  vector<string> reads;
  for (int c = 0; c < s.size() / 10000 * 40; c++) {
    string sc = Screw(s.substr(rand()%(s.size()-10000), 10000));
    reads.push_back(sc);
    vector<int> read_cluster;
    for (int i = 0; i + out < sc.size(); i+=100) {
      parts.push_back(DNASeq(sc.substr(i, out)));
      origin.push_back(c);
      origin_order.push_back(i / 100);
      read_cluster.push_back(-1);
      if (i > 0) {
        next.push_back(parts.size()-2);
      } else {
        next.push_back(-1);
      } 
    }
    read_clusters.push_back(read_cluster);
  }
  printf("parts size %d\n", parts.size());

  vector<int> order, cluster_id, rev(parts.size());
  vector<vector<unsigned int>> cluster_kmers(parts.size());
  vector<vector<unsigned int>> cluster_kmers2(parts.size());
  for (int i = 0; i < parts.size(); i++) {
    order.push_back(i); cluster_id.push_back(i);
  }

  random_shuffle(order.begin(), order.end());
  vector<unsigned int> cur_kmers, cur_kmers2;
  vector<bool> index_filter(1 << 28);
  vector<vector<int>*> index_parts;
  priority_queue<pair<int, pair<int,int>>, vector<pair<int, pair<int,int>>>, greater<pair<int,pair<int, int>>>> fr;

  for (int i = 0; i < order.size(); i++) {
    if (i % 100 == 0) {
      printf("i %d\r", i);
      fflush(stdout);
    }
    int el = order[i];
    rev[el] = i;
    DNASeq& ss = parts[el];
    int min_ord = i;
    cur_kmers.clear();
//    cur_kmers2.clear();
    PrepareKmers(ss, 14, cur_kmers);
//    PrepareKmers(ss, 6, cur_kmers2);
    index_parts.clear();
    while (!fr.empty()) fr.pop();
    for (auto &km: cur_kmers) {
      if (index_filter[km] == false) continue;
      index_parts.push_back(&index[km]);
    }
    for (int j = 0; j < index_parts.size(); j++) {
      fr.push(make_pair((*index_parts[j])[0], make_pair(j, 0)));
    }
    int last_cand = -1;
    int last_cand_count = 0;
    while (!fr.empty()) {
      auto x = fr.top();
      fr.pop();
      if ((*index_parts[x.second.first]).size() > x.second.second + 1) {
        fr.push(make_pair((*index_parts[x.second.first])[x.second.second+1],
                          make_pair(x.second.first, x.second.second+1)));
      }
      if (x.first != last_cand) {
        last_cand = x.first;
        last_cand_count = 0;
      }
      last_cand_count++;
      if (last_cand_count >= 10) {
        min_ord = last_cand;
        break;
      }
    }

/*    for (auto &km: cur_kmers) {
      if (index_filter[km] == false) continue;
      vector<int>& cands = index[km];
      for (int k = 0; k < cands.size(); k++) {
        if (cluster_id[cands[k]] != cands[k]) {
          continue;
        }
        int el_cand = order[cands[k]];
        if (origin[el] == origin[el_cand] && abs(origin_order[el] - origin_order[el_cand]) < 7) { 
          min_ord = min(min_ord, cands[k]);
          break;
        } else if (cands[k] < min_ord) {
          if (Check(cur_kmers, cluster_kmers[cands[k]], 14, 0.0025)) {
            if (true || Check(cur_kmers2, cluster_kmers2[cands[k]], 6, 0.1)) {
              min_ord = cands[k];
              break;
            } else {
            }
          }
        } else {
          break;
        }
      }
    }*/
    assert(cluster_id[min_ord] == min_ord);
    cluster_id[i] = min_ord;
    read_clusters[origin[el]][origin_order[el]] = cluster_id[i];
    if (min_ord == i) {
      cur_kmers2.clear();
      PrepareKmers(ss, 6, cur_kmers2);
      radix_sort2(cur_kmers.data(), cur_kmers.data() + cur_kmers.size());
      radix_sort2(cur_kmers2.data(), cur_kmers2.data() + cur_kmers2.size());
      cluster_kmers[i] = cur_kmers;
      cluster_kmers2[i] = cur_kmers2;
      for (auto &km: cur_kmers) {
        index[km].push_back(i);
        index_filter[km] = true;
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

  ProcessGraph(read_clusters, edges, cluster_sizes, cluster_kmers2, reads);
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
