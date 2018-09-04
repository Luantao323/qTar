#ifndef STP_H_
#define STP_H_

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <stdlib.h>
#include <map>
#include <vector>
#include <string>
#include <pthread.h>
#include "gfun.h"
#include "pssw/ssw_cpp.h"

using namespace std;

struct tarInf {
  string nam;
  int pos;
} tarStr;

struct rnaInf {
  string nam;
  string seq;
  int cmp;
} rnaStr;

struct thr {
  pthread_t thread;
  int id;
};

typedef pair<int, int> PINT;
int _pcmp (const PINT &x, const PINT &y) 
{
  return x.second > y.second;
}

string srna, target, output;
int thread = 1, length = 5, seed_num = 2, seed_beg = 2, seed_end = 8, mir_max = 5;
float score_aln = 18, score_fe = -10, seed_pen = 2;
bool help;
map<string, vector<tarInf> > tarPos;
//map<string, string> tarPos; // less memory, but longer running time.
map<int, vector<rnaInf> > rnaSeq;
map<string, string> tarSeq;
map<char, int> alnCigar;
map<string, float> pfe;
ofstream otfile;
pthread_mutex_t io_mutex;

#endif
