#include "dna_seq.h"
#include <iostream>
#include <deque>

using namespace std;

const int kKmerSize = 12;

void BuildIndex(const DNASeq& genome, vector<vector<int>>& index) {
  vector<pair<unsigned int, int>> mhs;
  genome.GetMinHashesWithPos(kKmerSize, 101, Hasher(), mhs);
  for (auto &e: mhs) {
    if (index.size() <= e.first) {
      index.resize(e.first+47);
    }
    assert(index.size() > e.first);
    if (index[e.first].size() > 0) {
      if (index[e.first].back() == e.second) {
        continue;
      }
    }
    index[e.first].push_back(e.second);
  }
}

inline void PushIfNotVisited(
    int dist, int cur_genome_pos, int cur_read_pos,
    int read_pos, int genome_pos, int iteration,
    deque<pair<int, pair<int, int>>>&fr, vector<vector<int>>& visited) {
  int gp = cur_genome_pos - genome_pos + read_pos + 20;
  if (visited[cur_read_pos + 1][gp] != iteration) {
    fr.push_back(make_pair(dist, make_pair(cur_genome_pos, cur_read_pos)));
    visited[cur_read_pos + 1][gp] = iteration;
  }
}

int ProcessHit(int genome_pos, int read_pos, const DNASeq& read, const DNASeq& genome) {
  static deque<pair<int, pair<int, int>>> fr;
  static int iteration = 0;
  iteration++;
  static vector<vector<int>> visited(read.size() + 47, vector<int>(read.size() + 47));
  assert(read.ExtractKmer(read_pos, kKmerSize) == genome.ExtractKmer(genome_pos, kKmerSize));
  // Forward
  int forward_errs = -1;
  fr.push_back(make_pair(0, make_pair(genome_pos + kKmerSize, read_pos + kKmerSize)));
  while (!fr.empty()) {
    pair<int, pair<int, int>> x = fr.front();
    fr.pop_front();
    if (x.first > 6) continue;
    if (x.second.second == read.size()) {
      forward_errs = x.first;
      fr.clear();
      break;
    }
    if (genome[x.second.first] == read[x.second.second]) {
      if (x.second.first + 1 < genome.size() || x.second.second + 1 == read.size()) {
        PushIfNotVisited(x.first, x.second.first + 1, x.second.second + 1,
                         read_pos, genome_pos, iteration, fr, visited);
      }
    } else {
      if (x.second.first + 1 < genome.size()) {
        PushIfNotVisited(x.first + 1, x.second.first + 1, x.second.second + 1,
                         read_pos, genome_pos, iteration, fr, visited);
        PushIfNotVisited(x.first + 1, x.second.first + 1, x.second.second,
                         read_pos, genome_pos, iteration, fr, visited);
      }
      PushIfNotVisited(x.first + 1, x.second.first, x.second.second + 1,
                       read_pos, genome_pos, iteration, fr, visited);
    }
  }
  if (forward_errs == -1) return -1;
  // Backward
  int backward_errs = -1;
  fr.push_back(make_pair(0, make_pair(genome_pos - 1, read_pos - 1)));
  while (!fr.empty()) {
    pair<int, pair<int, int>> x = fr.front();
    fr.pop_front();
    if (x.first > 6) continue;
    if (x.second.second == -1) {
      backward_errs = x.first;
      fr.clear();
      break;
    }
    if (genome[x.second.first] == read[x.second.second]) {
      if (x.second.first - 1 >= 0 || x.second.second - 1 == -1) {
        PushIfNotVisited(x.first, x.second.first - 1, x.second.second - 1,
                         read_pos, genome_pos, iteration, fr, visited);
      }
    } else {
      if (x.second.first - 1 >= 0) {
        PushIfNotVisited(x.first + 1, x.second.first - 1, x.second.second - 1,
                         read_pos, genome_pos, iteration, fr, visited);
        PushIfNotVisited(x.first + 1, x.second.first - 1, x.second.second,
                         read_pos, genome_pos, iteration, fr, visited);
      }
      PushIfNotVisited(x.first + 1, x.second.first, x.second.second - 1,
                       read_pos, genome_pos, iteration, fr, visited);
    }
  }
  if (backward_errs == -1) return -1;
  return backward_errs + forward_errs;
}

int main(int argc, char**argv) {
  vector<DNASeq> genomes;
  LoadFasta(argv[1], genomes);
  DNASeq genome = genomes[0];
  cout << "genome loaded " << genome.size() << endl;
  vector<vector<int>> index;
  BuildIndex(genome, index);
  cout << "index size " << index.size() << endl;
  int ss = 0;
  for (auto &e: index) ss += e.size();
  cout << "index entries " << ss << endl;

  FastqLoader reads_reader(argv[2]);
  DNASeq read;
  DNASeq read_rev;
  int hits = 0;
  int with_hit = 0;
  int reads = 0;
  while (reads_reader.Next(read, read_rev)) {
    reads++;
    bool hit = false;
    pair<unsigned int, int> mh = read.GetMinHashWithPos(kKmerSize, Hasher());
    if (mh.first < index.size()) {
      for (auto &genome_pos: index[mh.first]) {
        if (ProcessHit(genome_pos, mh.second, read, genome)) {
          hits += 1;
          hit = true;
        }
      }
    }
    pair<unsigned int, int> mhr = read_rev.GetMinHashWithPos(kKmerSize, Hasher());
    if (mhr.first < index.size()) {
      for (auto &genome_pos: index[mhr.first]) {
        if (ProcessHit(genome_pos, mhr.second, read_rev, genome)) {
          hits += 1;
          hit = true;
        }
      }
    }
    if (hit) with_hit++;
  }
  cout << hits << " " << with_hit << " " << reads << endl;
}
