#include <iostream>

using namespace std;

#include "kmc_api/kmc_file.h"

int main(int argc, char **argv) {
  CKMCFile kmer_db;
  uint32 klen, mode, min_counter, pref_len, sign_len, min_c, counter;
  uint64 tot_kmers, max_c;
  if (!kmer_db.OpenForListing(argv[1])) {
    cerr << "ERROR: cannot open " << argv[1] << endl;
    return 1;
  }
  kmer_db.Info(klen, mode, min_counter, pref_len, sign_len, min_c, max_c,
               tot_kmers);

  cout << tot_kmers << endl;

  return 0;
}
