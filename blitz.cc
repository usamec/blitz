#include <chrono>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <unordered_map>

using namespace std::chrono;
using namespace std;

const int kIndexKmer = 11;
char trans[256];

inline char ReverseBase(char a) {
  if (a == 'A') return 'T';
  if (a == 'C') return 'G';
  if (a == 'G') return 'C';
  if (a == 'T') return 'A';
}

inline string ReverseSeq(const string& x) {
  string ret;
  for (int i = x.length()-1; i >= 0; i--) {
    ret += ReverseBase(x[i]);
  }
  return ret;
}

const int nHashes = 10;

unsigned int hc[] = {
  0x0faaffaa,
  0xf0ffaaff,
  0x1e000000,
  0xe1444444,
  0x2dbb66bb,
  0xd2aaffaa,
  0x3cffaaff,
  0xc388dd88,
  0x4b222222,
  0xb4bbccbb
};


unsigned int HashKmer(unsigned int x, unsigned int cc) {
  return cc ^ (x << 7) ^ (x << 4);
}

void GetMinHashForSeq(const string& seq, vector<unsigned int>& hashes) {
  hashes.clear();
  hashes.resize(nHashes);
  unsigned int curhash = 0;
  for (int i = 0; i < kIndexKmer; i++) {
    curhash <<= 2;
    curhash += trans[seq[i]];
  }
  for (int j = 0; j < nHashes; j++) {
    hashes[j] = max(hashes[j], HashKmer(curhash, hc[j]));
  }
  for (int i = kIndexKmer; i < seq.length(); i++) {
    curhash <<= 2;
    curhash &= (1 << (2*kIndexKmer)) - 1;
    curhash += trans[seq[i]];
    for (int j = 0; j < nHashes; j++) {
      hashes[j] = max(hashes[j], HashKmer(curhash, hc[j]));
    }
  }  
}

void BuildIndex(
    const string& seq, int read_len,
    vector<unordered_map<unsigned int, vector<int>>>& index) {
  index.clear();
  index.resize(nHashes);
  vector<unsigned int> hashes;
  for (int i = 0; i < seq.length() - read_len + 1; i++) {
    string s = seq.substr(i, read_len);
    GetMinHashForSeq(s, hashes);
    for (int j = 0; j < nHashes; j++) {
      index[j][hashes[j]].push_back(i);
    }
  }
}
void BuildIndex2(
    const string& seq, int read_len,
    vector<unordered_map<unsigned int, vector<pair<int, int>>>>& index) {
  index.clear();
  index.resize(nHashes);
  vector<unsigned int> hashes;
  vector<unsigned int> last_hashes(nHashes);
  vector<int> int_start(nHashes);
  for (int i = 0; i < seq.length() - read_len + 1; i++) {
    string s = seq.substr(i, read_len);
    GetMinHashForSeq(s, hashes);
    for (int j = 0; j < nHashes; j++) {
      if (i == 0) {
        last_hashes[j] = hashes[j];
        int_start[j] = 0;
      } else {
        if (hashes[j] != last_hashes[j]) {
          index[j][last_hashes[j]].push_back(make_pair(int_start[j], i));
          last_hashes[j] = hashes[j];
          int_start[j] = i;
        }
      }
    }
  }
  for (int j = 0; j < nHashes; j++) {
    index[j][hashes[j]].push_back(make_pair(int_start[j], seq.length() - read_len + 1));
  }
}

void LoadGenome(char* filename, string& genome) {
  ifstream f(filename);
  string l;
  while (getline(f, l)) {
    if (l[0] == '>') continue;
    genome += l;
  }
}

void ShowElapsedTime(chrono::time_point<std::chrono::system_clock>& start) {
  chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
}

void Go(char *fn, vector<unordered_map<unsigned int, vector<int>>>& index, string& genome) {
  ifstream f(fn);
  string l, l2;
  int hits = 0;
  int with_hit = 0;
  int lines = 0;
  vector<unsigned int> hashes;
  vector<int> counts(genome.length());
  vector<bool> processed(genome.length());
  while (getline(f, l)) {
    getline(f, l2);
    getline(f, l);
    getline(f, l);
    string lr = ReverseSeq(l2);
    bool hit = false;
    GetMinHashForSeq(l2, hashes);
    for (int j = 0; j < nHashes; j++) {
      const vector<int>& hs = index[j][hashes[j]];
      hits += hs.size();
      for (int k = 0; k < hs.size(); k++) {
        counts[hs[k]]++;
        processed[hs[k]] = true;
      }
    }
    for (int j = 0; j < nHashes; j++) {
      const vector<int>& hs = index[j][hashes[j]];
      for (int k = 0; k < hs.size(); k++) {
        if (processed[hs[k]]) {
          if (counts[hs[k]] > 1) {
            hit = true;
          }
          processed[hs[k]] = false;
        }
        counts[hs[k]]--;
      }
    }

    GetMinHashForSeq(lr, hashes);
    for (int j = 0; j < nHashes; j++) {
      const vector<int>& hs = index[j][hashes[j]];
      hits += hs.size();
      for (int k = 0; k < hs.size(); k++) {
        counts[hs[k]]++;
        processed[hs[k]] = true;
      }
    }
    for (int j = 0; j < nHashes; j++) {
      const vector<int>& hs = index[j][hashes[j]];
      for (int k = 0; k < hs.size(); k++) {
        if (processed[hs[k]]) {
          if (counts[hs[k]] > 1) {
            hit = true;
          }
          processed[hs[k]] = false;
        }
        counts[hs[k]]--;
      }
    }
    if (hit) {
      with_hit++;
    }
    lines++;
  }
  cout << "hits " << hits << " lines " << lines << endl;
  cout << "with hit " << with_hit << endl;
}

void Go2(char *fn, vector<unordered_map<unsigned int, vector<pair<int, int>>>>& index,
         string& genome, char* output_file_name) {
  ifstream f(fn);
  ofstream of(output_file_name);
  string l, l2, l3;
  int hits = 0;
  int with_hit = 0;
  int lines = 0;
  vector<unsigned int> hashes;
  vector<int> counts(genome.length());
  vector<bool> processed(genome.length());
  vector<pair<int, int>> events;
  int depth;
  int req_depth = 3;
  while (getline(f, l)) {
    getline(f, l2);
    getline(f, l3);
    getline(f, l3);
    string lr = ReverseSeq(l2);
    bool hit = false;
    GetMinHashForSeq(l2, hashes);
    events.clear();
    int last_hit = -47;
    for (int j = 0; j < nHashes; j++) {
      if (index[j].count(hashes[j]) == 0) continue;
      const vector<pair<int, int>>& hs = index[j][hashes[j]];
      for (int k = 0; k < hs.size(); k++) {
        events.push_back(make_pair(hs[k].first, 1));
        events.push_back(make_pair(hs[k].second, -1));
      }
    }
    sort(events.begin(), events.end());
    depth = 0;
    for (int i = 0; i < events.size(); i++) {
      depth += events[i].second;
      if (depth >= req_depth) {
        hit = true;
        hits++;
        if (events[i].first - last_hit > 101) {
          last_hit = i;
          of << l << " " << events[i].first << "\n";
        }
      }
    }

    GetMinHashForSeq(lr, hashes);
    events.clear();
    last_hit = -47;
    for (int j = 0; j < nHashes; j++) {
      if (index[j].count(hashes[j]) == 0) continue;
      const vector<pair<int, int>>& hs = index[j][hashes[j]];
      for (int k = 0; k < hs.size(); k++) {
        events.push_back(make_pair(hs[k].first, 1));
        events.push_back(make_pair(hs[k].second, -1));
      }
    }
    sort(events.begin(), events.end());
    depth = 0;
    for (int i = 0; i < events.size(); i++) {
      depth += events[i].second;
      if (depth >= req_depth) {
        hit = true;
        hits++;
        if (events[i].first - last_hit > 101) {
          last_hit = i;
          of << l << " " << events[i].first << " R\n";
        }
      }
    }
    if (hit) {
      with_hit++;
    }
    lines++;
  }
  cout << "hits " << hits << " lines " << lines << endl;
  cout << "with hit " << with_hit << endl;
}


int main(int argc, char** argv) {
  trans['A'] = 1;
  trans['T'] = 2;
  trans['C'] = 3;
  trans['G'] = 0;
  chrono::time_point<std::chrono::system_clock> start_time = high_resolution_clock::now();

  string genome;
  LoadGenome(argv[1], genome);
  cout << "genome length " << genome.length() << endl;
  ShowElapsedTime(start_time);

  vector<unordered_map<unsigned int, vector<pair<int, int>>>> index; 
//  vector<unordered_map<unsigned int, vector<int>>> index; 
  BuildIndex2(genome, 101, index);
  ShowElapsedTime(start_time);

  Go2(argv[2], index, genome, argv[3]);
  cout << "total ";
  ShowElapsedTime(start_time);
}
