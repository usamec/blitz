#ifndef DNA_SEQ_H__
#define DNA_SEQ_H__

#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <fstream>

using namespace std;

class Hasher {
  unsigned int h;
 public:
  Hasher() : h(0x9b7cd3) {
  }
  Hasher(unsigned int h) : h(h) {}
  unsigned int operator()(unsigned int x) {
    unsigned int hh = x ^ h;
    hh ^= hh >> 16;
    hh *= 0x85ebca6b;
    hh ^= hh >> 13;
    hh *= 0xc2b2ae35;
    hh ^= hh >> 16;
    return hh;
  }
};

class DNASeq {
  vector<unsigned int> storage_;
  int size_;
  static char trans[256];
 public:
  static void InitTrans() {
    trans['A'] = trans['a'] = 0;
    trans['C'] = trans['c'] = 1;
    trans['G'] = trans['g'] = 2;
    trans['T'] = trans['t'] = 3;
  }
  DNASeq() {}
  DNASeq(string x) {
    size_ = x.size();
    storage_.resize((size_ + 15) / 16 + 1);
    for (int i = 0; i < x.size(); i++) {
      storage_[i/16] |= (trans[x[i]] << (i % 16 * 2));
    }
  }

  char operator[](int i) const {
    return (storage_[i/16] >> (i % 16 * 2)) & 3;
  }

  static char BaseToStr(char b) {
    char t[] = "ACGT";
    return t[b];
  }

  static string KmerToStr(unsigned int kmer, int k) {
    string ret;
    for (int i = 0; i < k; i++) {
      ret += BaseToStr(kmer & 3);
      kmer >>= 2;
    }
    return ret;
  }

  // Works for k up to 16
  unsigned int ExtractKmer(int i, int k) const {
    assert(i+k <= size_);
    unsigned long join = ((unsigned long) storage_[i/16 + 1] << 32) | storage_[i/16];
    return (join >> (i % 16 * 2)) & ((1ll << (2*k)) - 1);
  }

  int size() const {
    return size_;
  }

  unsigned int GetMinHash(int kmer_size, Hasher hasher) const {
    unsigned int minhash = hasher(ExtractKmer(0, kmer_size));
    for (int i = 1; i < size_ - kmer_size + 1; i++) {
      minhash = min(minhash, hasher(ExtractKmer(i, kmer_size)));
    }
    return minhash;
  }
  
  unsigned int KillMiddle(unsigned int x) const {
    unsigned int masklow = 0xffff;
    unsigned int maskhigh = 0xffffffff << 18;
    return (x & masklow) | ((x & maskhigh) >> 2);
  }

  unsigned int GetMinHash2(int kmer_size, Hasher hasher) const {
    unsigned int minhash = hasher(ExtractKmer(0, kmer_size));
    for (int i = 1; i < size_ - kmer_size + 1; i++) {
      minhash = min(minhash, hasher(ExtractKmer(i, kmer_size)));
    }
    for (int i = 1; i < size_ - kmer_size; i++) {
      minhash = min(minhash, hasher(KillMiddle(ExtractKmer(i, kmer_size+1))));
    }
    return minhash;
  }

  pair<unsigned int, int> GetMinHashWithPos(int kmer_size, Hasher hasher) const {
    unsigned int minhash = hasher(ExtractKmer(0, kmer_size));
    int pos = 0;
    for (int i = 1; i < size_ - kmer_size + 1; i++) {
      unsigned int h = hasher(ExtractKmer(i, kmer_size));
      if (h < minhash) {
        minhash = h;
        pos = i;
      }
    }
    return make_pair(minhash, pos);
  }

  void GetMinHashesSlow(int kmer_size, int minhash_size, Hasher hasher,
                        vector<unsigned int>& minhashes) const {
    for (int i = 0; i < size_ - minhash_size + 1; i++) {
      unsigned int minhash = hasher(ExtractKmer(i, kmer_size)); 
      for (int j = 1; j < minhash_size - kmer_size + 1; j++) {
        minhash = min(minhash, hasher(ExtractKmer(i+j, kmer_size)));
      }
      minhashes.push_back(minhash);
    }
  }

  void GetMinHashes(int kmer_size, int minhash_size, Hasher hasher,
                    vector<unsigned int>& minhashes) const {
    vector<unsigned int> block1;
    vector<unsigned int> block2;
    vector<unsigned int> suffix_min;
    vector<unsigned int> prefix_min;
    int jump = minhash_size - kmer_size + 1;
    for (int i = 0; i < jump; i++) {
      block1.push_back(hasher(ExtractKmer(i, kmer_size)));
    }
    unsigned first_mh = block1[0];
    for (int i = 1; i < jump; i++) {
      first_mh = min(first_mh, block1[i]);
    }
    minhashes.push_back(first_mh);
    int cur_pos = jump;
    while (cur_pos < size_) {
      block2.clear();
      for (int i = cur_pos; i < cur_pos + jump && i < size_ - kmer_size + 1; i++) {
        block2.push_back(hasher(ExtractKmer(i, kmer_size)));
      }
      suffix_min.clear();
      prefix_min.clear();
      suffix_min.push_back(0xffffffff);
      for (int i = jump - 1; i >= 0; i--) {
        suffix_min.push_back(min(suffix_min.back(), block1[i]));
      }

      prefix_min.push_back(block2[0]);
      for (int i = 1; i < block2.size(); i++) {
        prefix_min.push_back(min(prefix_min.back(), block2[i]));
      }

      for (int i = 0; i < block2.size(); i++) {
        minhashes.push_back(min(prefix_min[i], suffix_min[jump - i - 1]));
      }

      swap(block1, block2);
      cur_pos += jump;
    }
  }

  void GetMinHashesWithPos(int kmer_size, int minhash_size, Hasher hasher,
                    vector<pair<unsigned int, int>>& minhashes) const {
    vector<pair<unsigned int,int>> block1;
    vector<pair<unsigned int,int>> block2;
    vector<pair<unsigned int,int>> suffix_min;
    vector<pair<unsigned int,int>> prefix_min;
    int jump = minhash_size - kmer_size + 1;
    for (int i = 0; i < jump; i++) {
      block1.push_back(make_pair(hasher(ExtractKmer(i, kmer_size)), i));
    }
    pair<unsigned int,int> first_mh = block1[0];
    for (int i = 1; i < jump; i++) {
      first_mh = min(first_mh, block1[i]);
    }
    minhashes.push_back(first_mh);
    int cur_pos = jump;
    while (cur_pos < size_) {
      block2.clear();
      for (int i = cur_pos; i < cur_pos + jump && i < size_ - kmer_size + 1; i++) {
        block2.push_back(make_pair(hasher(ExtractKmer(i, kmer_size)), i));
      }
      suffix_min.clear();
      prefix_min.clear();
      suffix_min.push_back(make_pair(0xffffffff, size_+47));
      for (int i = jump - 1; i >= 0; i--) {
        suffix_min.push_back(min(suffix_min.back(), block1[i]));
      }

      prefix_min.push_back(block2[0]);
      for (int i = 1; i < block2.size(); i++) {
        prefix_min.push_back(min(prefix_min.back(), block2[i]));
      }

      for (int i = 0; i < block2.size(); i++) {
        minhashes.push_back(min(prefix_min[i], suffix_min[jump - i - 1]));
      }

      swap(block1, block2);
      cur_pos += jump;
    }
  }
};

void LoadFasta(char* filename, vector<DNASeq>& seqs);
void LoadFastq(char* filename, vector<DNASeq>& seqs);

class FastqLoader {
  ifstream f;
  string s, s2, sr;
  char rev[256];
  void ReverseSeq(const string& x, string& ret) {
    ret.clear();
    for (int i = x.length()-1; i >= 0; i--) {
      ret += rev[x[i]];
    }
  }
 public:
  FastqLoader(char* filename) : f(filename) {
    rev['A'] = 'T'; rev['T'] = 'A';
    rev['C'] = 'G'; rev['G'] = 'C';
  }

  bool Next(DNASeq& seq, DNASeq& seq_rev, string& name) {
    if (!getline(f, name)) return false;
    getline(f, s);
    getline(f, s2);
    getline(f, s2);
    seq = s;
    ReverseSeq(s, sr);
    seq_rev = sr;
    return true;
  }
};

#endif
