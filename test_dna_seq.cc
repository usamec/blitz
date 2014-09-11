#include "dna_seq.h"
#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <chrono>
using namespace std::chrono;

void ShowElapsedTime(chrono::time_point<std::chrono::system_clock>& start) {
  chrono::time_point<std::chrono::system_clock> end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  cout << "elapsed time: " << elapsed_seconds.count() << "s\n";
}
int main() {
  char a[] = "ACGT";
  string ss;
  srand(time(NULL));
  for (int i = 0; i < 20000000; i++) {
    ss += a[rand()%4];
  }
  DNASeq s(ss);
  for (int i = 0; i < ss.size(); i++) {
    assert(ss[i] == DNASeq::BaseToStr(s[i]));
  }

  for (int i = 0; i < ss.size()-13 && i < 50; i++) {
    string s1 = ss.substr(i, 13);
    string s2 = DNASeq::KmerToStr(s.ExtractKmer(i, 13), 13);
    assert(s1 == s2);
  }

  chrono::time_point<std::chrono::system_clock> start_time = high_resolution_clock::now();
  vector<unsigned int> mh, mh2;
  s.GetMinHashesSlow(12, 100, Hasher(), mh);
  cout << mh.size() << endl;
  ShowElapsedTime(start_time);
  s.GetMinHashes(12, 100, Hasher(), mh2);
  cout << mh2.size() << endl;
  ShowElapsedTime(start_time);
  assert(mh.size() == mh2.size());
  for (int i = 0; i < mh.size(); i++) {
    if (mh[i] != mh2[i]) {
      printf("fuck %d %X %X\n", i, mh[i], mh2[i]);
    }
  }
}
