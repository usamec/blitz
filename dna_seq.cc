#include "dna_seq.h"
#include <fstream>
#include <iostream>

using namespace std;

void LoadFasta(char* filename, vector<DNASeq>& seqs) {
  ifstream f(filename);
  string s;
  string buf;
  while (getline(f, s)) {
    if (s[0] == '>') {
      if (buf.size() > 0) {
        seqs.push_back(DNASeq(buf));
        buf = "";
      }
    } else {
      buf += s;
    }
  }
  if (buf.size() > 0) {
    seqs.push_back(DNASeq(buf));
  }
}

void LoadFastq(char* filename, vector<DNASeq>& seqs) {
  ifstream f(filename);
  string s, s2;
  while (getline(f, s2)) {
    getline(f, s);
    getline(f, s2);
    getline(f, s2);
    seqs.push_back(DNASeq(s));
  }
}
