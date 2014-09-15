#ifndef GRAPH_MANIPULATOR_H__
#define GRAPH_MANIPULATOR_H__

#include <cstdio>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <vector>

const int kCutoff = 2;

class GraphManipulator {
 public:
  GraphManipulator(vector<vector<int>>& read_clusters, 
                   map<pair<int, int>, int>& edges,
                   unordered_map<int, int>& cluster_sizes) : 
      read_clusters(read_clusters), cluster_sizes(cluster_sizes), edge_cov(edges) {
    BuildInitialGraph(edges);
    for (auto &e: read_clusters) {
      vector<int> rc;
      int last = -1;
      for (auto &e2: e) {
        if (e2 != last) rc.push_back(e2);
        last = e2;
      }
      read_clusters_compact.push_back(rc);
    }
  }
  
  void BuildInitialGraph(map<pair<int, int>, int>& edges) {
    unordered_map<int, int> trans;
    for (auto &e: edges) {
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

  bool RemoveTips() {
    bool changed = false;
    for (int i = 0; i < g.size(); i++) {
      if (g[i].size() + gr[i].size() == 1) {
        int total_cov = 0;
        for (int j = 0; j < in_node_seq[i].size(); j++) {
          total_cov += cluster_sizes[in_node_seq[i][j]];
        }
        if (total_cov < 20) {
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
    }
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
//    printf("rc %d %d %d\n", start, read_cluster[start], node_seq.size());
    bool match = 0;
    for (int i = 0; i < node_seq.size(); i++) {
      if (start - 1 - i < 0) return -1;
//      printf("  test %d %d %d\n", i, node_seq[i], read_cluster[start+1+i]);
      if (node_seq[i] != read_cluster[start-1-i]) return -1;
    }
    if (start - 1 - node_seq.size() < 0) return -1;
    return read_cluster[start-1-node_seq.size()];
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
    for (int i = 0; i < g.size(); i++) {
      if (in_node_seq[i].size() == 0) continue;
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

      }
    }

    return changed;
  }

  void Process() {
    printf("g size %d\n", g.size());

    bool changed = true;
    while (changed) {
      changed = false;

      changed |= RemoveTips();
      changed |= RemoveBubbles();
      changed |= Concat();
//      changed |= ThreadReads();
      changed |= ThreadReads2();

      printf("g size %d\n", GetSize());
    }
  }

  int GetEdgeCov(int from, int to) {
    int cluster_from = in_node_seq[from].back();
    int cluster_to = in_node_seq[to][0];
    return edge_cov[make_pair(cluster_from, cluster_to)];
  }

  void OutputGraph() {
    FILE *f = fopen("wtf.dot", "w");
    fprintf(f, "digraph {\n");
    for (int i = 0; i < in_node_seq.size(); i++) {
      if (g[i].size() == 0 && gr[i].size() == 0) continue;
      fprintf(f, "%d [label=\"", i);
      for (int j = 0; j < in_node_seq[i].size(); j++) {
        fprintf(f, "%d (%d)%c", in_node_seq[i][j], cluster_sizes[in_node_seq[i][j]],
                j+1 == in_node_seq[i].size() ? '"' : ',');
      }
      fprintf(f, "]\n");
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
};

#endif 
