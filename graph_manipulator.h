#ifndef GRAPH_MANIPULATOR_H__
#define GRAPH_MANIPULATOR_H__

#include <cstdio>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <vector>
#include <queue>

const int kCutoff = 5;

class GraphManipulator {
 public:
  GraphManipulator(vector<vector<int>>& read_clusters, 
                   map<pair<int, int>, int>& edges,
                   unordered_map<int, int>& cluster_sizes,
                   vector<vector<unsigned int>>& cluster_kmers) : 
      read_clusters(read_clusters), cluster_sizes(cluster_sizes), edge_cov(edges),
      cluster_kmers(cluster_kmers) {
    for (auto &e: read_clusters) {
      vector<int> rc;
      int last = -1;
      for (auto &e2: e) {
        if (e2 != last) rc.push_back(e2);
        last = e2;
      }
      read_clusters_compact.push_back(rc);
    }
    BuildInitialGraph();
  }
  
  void BuildInitialGraph() {
    unordered_map<int, int> trans;
    RecalculateEdgeCov();
    for (auto &e: edge_cov) {
      if (e.second <= kCutoff) continue;
      if (trans.count(e.first.first) == 0) {
        int id = trans.size();
        g.push_back(unordered_set<int>());
        gr.push_back(unordered_set<int>());
        in_node_seq.push_back(vector<int>({e.first.first}));
        trans[e.first.first] = id;
      }
      if (trans.count(e.first.second) == 0) {
        int id = trans.size();
        g.push_back(unordered_set<int>());
        gr.push_back(unordered_set<int>());
        in_node_seq.push_back(vector<int>({e.first.second}));
        trans[e.first.second] = id;
      }
      g[trans[e.first.first]].insert(trans[e.first.second]);
      gr[trans[e.first.second]].insert(trans[e.first.first]);
    }
  }

  bool IsCoveredElseAndHasSmallCov(int i) {
    if (g[i].size() > 0) {
      int parent = *g[i].begin();
      int end_v = in_node_seq[parent].back();
      int total_cov = 0;
      for (auto &e: g[parent]) {
        total_cov += edge_cov[make_pair(end_v, in_node_seq[e][0])];
      }
      int our_cov = edge_cov[make_pair(end_v, in_node_seq[i][0])];
      if (our_cov * 3 > total_cov) return false;     
    }
    unordered_set<int> missing(in_node_seq[i].begin(), in_node_seq[i].end());
    for (int j = 0; j < g.size(); j++) {
      if (i == j) continue;
      for (auto &e: in_node_seq[j]) {
        missing.erase(e);
      }
    }
    return missing.empty();
  }

  bool RemoveTips() {
    RecalculateEdgeCov();
    bool changed = false;
    for (int i = 0; i < g.size(); i++) {
      if (g[i].size() + gr[i].size() == 1) {
        int total_cov = 0;
        for (int j = 0; j < in_node_seq[i].size(); j++) {
          total_cov += cluster_sizes[in_node_seq[i][j]];
        }
        if (total_cov < 60 || IsCoveredElseAndHasSmallCov(i)) {
          changed = true;
          for (auto &next: g[i]) {
            gr[next].erase(i);
          }
          for (auto &prev: gr[i]) {
            g[prev].erase(i);
          }
          g[i].clear();
          gr[i].clear();
          in_node_seq[i].clear();
        }
      }
    }
    return changed;
  }

  void RemoveFromReads(int cluster_id) {
    for (auto &r: read_clusters) {
      vector<int> ok;
      for (int i = 0; i < r.size(); i++) {
        if (r[i] != cluster_id) ok.push_back(r[i]);
      }
      r = ok;
    }
  
    for (auto &r: read_clusters_compact) {
      vector<int> ok;
      for (int i = 0; i < r.size(); i++) {
        if (r[i] != cluster_id) ok.push_back(r[i]);
      }
      r = ok;
    }
    for (auto &r: in_node_seq) {
      vector<int> ok;
      for (int i = 0; i < r.size(); i++) {
        if (r[i] != cluster_id) ok.push_back(r[i]);
      }
      r = ok;
    }
  }

  void EnumerateUpToDist(int i, int dist, vector<int>& out) {
    queue<int> fr;
    unordered_map<int, int> d;
    d[i] = 0;
    fr.push(i);
    while (!fr.empty()) {
      int x = fr.front(); fr.pop();
      if (d[x] >= 2) out.push_back(x);
      if (d[x] >= dist) continue;
      for (auto &n: g[x]) {
        if (d.count(n) != 0) continue;
        d[n] = d[x] + 1;
        fr.push(n);
      }
    }
  }

  void EnumerateUpToDistBlock(int i, int dist, unordered_set<int>& out, int block,
                              vector<unordered_set<int>>& gg) {
    queue<int> fr;
    unordered_map<int, int> d;
    d[i] = 0;
    fr.push(i);
    while (!fr.empty()) {
      int x = fr.front(); fr.pop();
      if (x != i && x != block) out.insert(x);
      if (x == block) continue;
      if (d[x] >= dist) continue;
      for (auto &n: gg[x]) {
        if (d.count(n) != 0) continue;
        d[n] = d[x] + 1;
        fr.push(n);
      }
    }
  }

  bool IsClosed(int i, int c, vector<int>& closure) {
    unordered_set<int> cf, cb;
    EnumerateUpToDistBlock(i, 5, cf, c, g);
    EnumerateUpToDistBlock(c, 5, cb, i, gr);
    if (cf != cb) return false;
    for (auto &e: cf) {
      if (in_node_seq[e].size() > 1) return false;
    }
    closure.clear();
    closure.insert(closure.end(), cf.begin(), cf.end());
    return true;
  }

  bool RemoveBubbles() {
    bool changed = false;
    for (int i = 0; i < g.size(); i++) {
      if (g[i].size() == 1 && gr[i].size() == 1 && in_node_seq[i].size() == 1) {
        int prev = *gr[i].begin();
        int next = *g[i].begin();
        if (prev == next) continue;
        if (g[prev].count(next)) {
          g[prev].erase(i);
          gr[next].erase(i);
          g[i].clear();
          gr[i].clear();
          RemoveFromReads(in_node_seq[i][0]);
          in_node_seq[i].clear();
          changed = true;
        }
      }
    }
    return changed;
  }

  bool RemoveBubbles2() {
    bool changed = false;
    for (int i = 0; i < g.size(); i++) {
      vector<int> cands;
      EnumerateUpToDist(i, 3, cands);
      for (auto &c: cands) {
        vector<int> closure;
        if (!IsClosed(i, c, closure)) continue;
        if (closure.size() < 2) continue;
        changed = true;
        g[i].insert(c);
        gr[c].insert(i);
        for (auto &e: closure) {
          g[i].erase(e);
          gr[c].erase(e);
          g[e].clear();
          gr[e].clear();
          assert(in_node_seq[e].size() <= 1);
          for (auto &s: in_node_seq[e]) {
            RemoveFromReads(s);     
          }
        }
      }
    }

    return changed;
  }

  bool Concat() {
    bool changed = false;
    for (int i = 0; i < g.size(); i++) {
      if (g[i].size() != 1) continue;
      int next = *g[i].begin();
      if (gr[next].size() != 1) continue;
      if (g[next].count(i) != 0) continue;
      g[i] = g[next];
      for (auto &e: g[next]) {
        gr[e].erase(next);
        gr[e].insert(i);
      }
      in_node_seq[i].insert(
          in_node_seq[i].end(), in_node_seq[next].begin(), in_node_seq[next].end());
      in_node_seq[next].clear();
      g[next].clear();
      gr[next].clear();
      changed = true;
    }
    return changed;
  }

  int GetSize() {
    int s = 0;
    for (int i = 0; i < g.size(); i++) {
      if (in_node_seq[i].size() > 0) s++;
    }
    return s;
  }
  
  void GetReadsWithCluster(int cluster, vector<pair<int, int>>& reads) {
    for (int i = 0; i < read_clusters.size(); i++) {
      for (int j = 0; j < read_clusters[i].size(); j++) {
        if (read_clusters[i][j] == cluster) {
          reads.push_back(make_pair(i, j));
          break;
        }
      }
    }
  }


  bool ThreadReads() {
    bool changed = false;

    for (int i = 0; i < g.size(); i++) {
      if (g[i].size() >= 2 && gr[i].size() >= 2 && in_node_seq[i].size() == 1) {
        vector<pair<int, int>> our_reads;
        GetReadsWithCluster(in_node_seq[i][0], our_reads);
        unordered_map<int, unordered_map<int, int>> before_to_after_counts;
        unordered_map<int, unordered_map<int, int>> after_to_before_counts;
        for (auto &r: our_reads) {
          if (r.second == 0) continue;
          int after = read_clusters[r.first][r.second-1];
          int before = -1;
          for (int j = r.second; j < read_clusters[r.first].size(); j++) {
            if (read_clusters[r.first][j] != in_node_seq[i][0]) {
              before = read_clusters[r.first][j];
              break;
            }
          }
          if (before == -1) continue;
          before_to_after_counts[before][after]++;  
          after_to_before_counts[after][before]++;  
        }
        bool all_one = true;
        for (auto &e: after_to_before_counts) {
          vector<int> out;
          for (auto &e2: e.second) {
            if (e2.second <= kCutoff) out.push_back(e2.first);
          }
          for (auto &e2: out) {
            e.second.erase(e2);
          }
        }
        for (auto &e: before_to_after_counts) {
          vector<int> out;
          for (auto &e2: e.second) {
            if (e2.second <= kCutoff) out.push_back(e2.first);
          }
          for (auto &e2: out) {
            e.second.erase(e2);
          }
          if (e.second.size() == 1) {
            if (e.second.begin()->first == e.first) continue;
            if (after_to_before_counts[e.second.begin()->first].size() != 1) continue;
            vector<int> prevs, nexts;
            for (auto &prev: gr[i]) {
              if (g[i].count(prev)) continue;
              if (g[prev].size() > 1) continue;
              if (in_node_seq[prev].back() == e.first) {
                prevs.push_back(prev);
              }
            }
            for (auto &next: g[i]) {
              if (gr[i].count(next)) continue;
              if (in_node_seq[next][0] == e.second.begin()->first) {
                nexts.push_back(next);
              }
            }
            printf("split %d %d %d\n", in_node_seq[i][0], e.first, e.second.begin()->first);
            if (prevs.size() == 0 || nexts.size() == 0)
              continue;
              
            for (auto &prev: prevs) {
              changed = true;
              g[prev].erase(i);
              gr[i].erase(prev);
              in_node_seq[prev].push_back(in_node_seq[i][0]);
              for (auto &next: nexts) {
                g[prev].insert(next);
                gr[next].erase(i);
                gr[next].insert(prev);
                g[i].erase(next); // TODO: check other direction
              }
            }
          }
        }
        if (g[i].size() == 0 && gr[i].size() == 0) {
          in_node_seq[i].clear();
        }
      }
    }

    return changed;
  }

  int GetClusterAfterSeq(vector<int>& read_cluster, int start, vector<int>& node_seq) {
    bool match = 0;
    for (int i = 0; i < node_seq.size(); i++) {
      if (start - 1 - i < 0) return -1;
      if (node_seq[i] != read_cluster[start-1-i]) return -1;
    }
    if (start - 1 - node_seq.size() < 0) return -1;
    return read_cluster[start-1-node_seq.size()];
  }

  int GetClusterBeforeSeq(vector<int>& read_cluster, int start, vector<int>& node_seq) {
    bool match = 0;
    for (int i = 0; i < node_seq.size(); i++) {
      if (start + 1 + i >= read_cluster.size()) return -1;
      if (node_seq[node_seq.size() - 1 - i] != read_cluster[start+1+i]) return -1;
    }
    if (start + 1 + node_seq.size() >= read_cluster.size()) return -1;
    return read_cluster[start+1+node_seq.size()];
  }

  void GetReadsWithCluster2(int cluster, vector<pair<int, int>>& reads) {
    for (int i = 0; i < read_clusters_compact.size(); i++) {
      for (int j = 0; j < read_clusters_compact[i].size(); j++) {
        if (read_clusters_compact[i][j] == cluster) {
          reads.push_back(make_pair(i, j));
          break;
        }
      }
    }
  }

  bool ThreadReads2() {
    bool changed = false;
    vector<vector<int>> forward;
    vector<vector<int>> backward;

    for (int i = 0; i < g.size(); i++) {
      if (in_node_seq[i].size() == 0) continue;
      if (g[i].size() != 1) continue;
      vector<pair<int, int>> our_reads;
      GetReadsWithCluster2(in_node_seq[i].back(), our_reads);
      for (auto &e: g[i]) {
        unordered_map<int, int> next_clusters;
        for (int j = 0; j < our_reads.size(); j++) {
          int next_cluster = GetClusterAfterSeq(
              read_clusters_compact[our_reads[j].first], our_reads[j].second, in_node_seq[e]);
          if (next_cluster == -1) continue;
          next_clusters[next_cluster]++;
        }
        vector<int> out;
        for (auto &e2: next_clusters) {
          if (e2.second <= kCutoff) out.push_back(e2.first);
        }
        for (auto &e2: out) next_clusters.erase(e2);
        if (next_clusters.size() != 1) continue;
        printf("go %d %d %d\n", in_node_seq[i].back(), in_node_seq[e].back(),
            next_clusters.begin()->first);
        int next_cluster = next_clusters.begin()->first;
        int next_cluster_node = -1;
        for (auto &e2: g[e]) {
          if (in_node_seq[e2][0] == next_cluster) {
            if (next_cluster_node == -1) next_cluster_node = e2;
            else {
              printf("fail\n");
              next_cluster_node = -2;
            }
          }
        }
        if (next_cluster_node >= 0) {
          forward.push_back(vector<int>({i, e, next_cluster_node}));
        }
      }
    }

    for (int i = 0; i < g.size(); i++) {
      if (in_node_seq[i].size() == 0) continue;
//      if (gr[i].size() != 1) continue;
      vector<pair<int, int>> our_reads;
      GetReadsWithCluster2(in_node_seq[i][0], our_reads);
      for (auto &e: gr[i]) {
        unordered_map<int, int> next_clusters;
        for (int j = 0; j < our_reads.size(); j++) {
          int next_cluster = GetClusterBeforeSeq(
              read_clusters_compact[our_reads[j].first], our_reads[j].second, in_node_seq[e]);
          if (next_cluster == -1) continue;
          next_clusters[next_cluster]++;
        }
        vector<int> out;
        for (auto &e2: next_clusters) {
          if (e2.second <= kCutoff) out.push_back(e2.first);
        }
        for (auto &e2: out) next_clusters.erase(e2);
        if (next_clusters.size() != 1) continue;
        printf("gob %d %d %d\n", in_node_seq[i][0], in_node_seq[e].back(),
            next_clusters.begin()->first);
        int next_cluster = next_clusters.begin()->first;
        int next_cluster_node = -1;
        for (auto &e2: gr[e]) {
          if (in_node_seq[e2].back() == next_cluster) {
            if (next_cluster_node == -1) next_cluster_node = e2;
            else { 
              printf("fail\n");
              next_cluster_node = -2;
            }
          }
        }
        printf("ncn %d\n", next_cluster_node);
        if (next_cluster_node >= 0) {
          backward.push_back(vector<int>({i, e, next_cluster_node}));
        }
      }
    }

    for (auto &f: forward) {
      for (auto &b: backward) {
        if (b[0] == f[2] && b[1] == f[1] && b[2] == f[0]) {
          g[f[0]].erase(f[1]);
          in_node_seq[f[0]].insert(in_node_seq[f[0]].end(),
              in_node_seq[f[1]].begin(), in_node_seq[f[1]].end());
          g[f[0]].insert(f[2]);
          gr[f[1]].erase(f[0]);
          g[f[1]].erase(f[2]);
          gr[f[2]].erase(f[1]);
          gr[f[2]].insert(f[0]);
          if (g[f[1]].size() == 0 && gr[f[1]].size() == 0) {
            in_node_seq[f[1]].clear();
          }
          return true;
        }
      }
    }

    return changed;
  }

  bool CanMerge(const vector<unsigned int>& ka, const vector<unsigned int>& kb) {
    static vector<unsigned int> inter;
    inter.clear();
    inter.resize(ka.size() + kb.size());
    auto it = set_intersection(ka.begin(), ka.end(), kb.begin(), kb.end(), inter.begin());
    int inter_size = distance(inter.begin(), it);
    double jaccard = 1.* inter_size / (ka.size() + kb.size() - inter_size);

    printf("%lf\n", jaccard);
    return jaccard >= 0.001;
  }

  bool MergeNeighbors() {
    bool changed = false;

    printf("testing merge\n");
    for (int i = 0; i < g.size(); i++) {
      for (auto &e: g[i]) {
        for (auto &e2: g[i]) {
          if (e == e2) continue;
          if (e > e2) continue;
          if (in_node_seq[e].size() > 1 || in_node_seq[e2].size() > 1) continue;
          printf("test %d %d ", in_node_seq[e][0], in_node_seq[e2][0]); 
          if (CanMerge(cluster_kmers[in_node_seq[e][0]], cluster_kmers[in_node_seq[e2][0]])) {
            printf("Merge %d %d\n", in_node_seq[e][0], in_node_seq[e2][0]);
          }
        }
      }
    }

    return changed;
  }

  void RemoveSmall() {
    for (int i = 0; i < g.size(); i++) {
      if (in_node_seq[i].size() != 1) continue;
      printf("remove %d\n", in_node_seq[i][0]);
      for (auto &prev: gr[i]) {
        for (auto &next: g[i]) {
          g[prev].insert(next);
          gr[next].insert(prev);
        }
        g[prev].erase(i);
      }
      for (auto &next: g[i]) {
        gr[next].erase(i);
      }
      RemoveFromReads(in_node_seq[i][0]);
      in_node_seq[i].clear();
      g[i].clear();
      gr[i].clear();
    }
    RecalculateEdgeCov();

    for (int i = 0; i < g.size(); i++) {
      int scov = 0;
      for (auto &e: g[i]) {
        scov += GetEdgeCov(i, e);
      }
      vector<int> out;
      for (auto &e: g[i]) {
        if (GetEdgeCov(i, e) == 0 || GetEdgeCov(i, e) < scov / 4) {
          out.push_back(e);
        }
      }
      for (auto &e: out) {
        g[i].erase(e);
        gr[e].erase(i);
      }
    }
    bool changed = true;
    while (changed) {
      changed = false;
      changed |= Concat();
    }
  }

  void ThreadBig() {
    RemoveSmall();
  }

  void Process() {
    printf("g size %d\n", g.size());

    bool changed = true;
    int id = 1;
    char buf[20];
    while (changed) {
      changed = false;

      changed |= RemoveTips();
      sprintf(buf, "wtf%dt.dot", id);
      OutputGraph(buf);
      changed |= RemoveBubbles();
      sprintf(buf, "wtf%db.dot", id);
      OutputGraph(buf);
      changed |= Concat();
      sprintf(buf, "wtf%dc.dot", id);
      OutputGraph(buf);
//      changed |= ThreadReads();
      changed |= ThreadReads2();
      sprintf(buf, "wtf%dth.dot", id);
      OutputGraph(buf);
      if (!changed) {
        changed |= RemoveBubbles2();
        sprintf(buf, "wtf%db2.dot", id);
        OutputGraph(buf);
      }

      id++;
      printf("g size %d\n", GetSize());
    }
    if (GetSize() > 1) {
      ThreadBig();
    }
    printf("done %d\n", GetSize());
  }

  int GetEdgeCov(int from, int to) {
    int cluster_from = in_node_seq[from].back();
    int cluster_to = in_node_seq[to][0];
    return edge_cov[make_pair(cluster_from, cluster_to)];
  }

  void RecalculateEdgeCov() {
    edge_cov.clear();
    for (auto &e: read_clusters_compact) {
      for (int i = 1; i < e.size(); i++) {
        edge_cov[make_pair(e[i], e[i-1])]++;
      }
    }
  }

  void OutputGraph(char* filename, bool all=true) {
    RecalculateEdgeCov();
    FILE *f = fopen(filename, "w");
    fprintf(f, "digraph {\n");
    for (int i = 0; i < in_node_seq.size(); i++) {
      if (in_node_seq[i].size() == 0) continue;
      fprintf(f, "%d [label=\"", i);
      for (int j = 0; j < in_node_seq[i].size(); j++) {
        fprintf(f, "%d (%d)%c", in_node_seq[i][j], cluster_sizes[in_node_seq[i][j]],
                j+1 == in_node_seq[i].size() ? ' ' : ',');
        if (!all) {
          j += max(0, (int)in_node_seq[i].size() - 2);
        }
      }
      if (!all) 
        fprintf(f, "s %d\"]\n", in_node_seq[i].size());
      else
        fprintf(f, "\"]\n");
    }
    for (int i = 0; i < g.size(); i++) {
      for (auto &e: g[i]) {
        fprintf(f, "%d -> %d [label=\"%d\"]\n", i, e, GetEdgeCov(i, e));
      }
    }
    fprintf(f, "}\n");
    fclose(f);
  }

  vector<vector<int>> read_clusters;
  vector<vector<int>> read_clusters_compact;
  unordered_map<int, int> cluster_sizes;
  vector<unordered_set<int>> g, gr;
  vector<vector<int>> in_node_seq;
  map<pair<int, int>, int> edge_cov;
  vector<vector<unsigned int>>& cluster_kmers;
};

#endif 
