#include "sRNATarPredictor.h"

void _getopt(int argc, char* argv[]);
void _usage();
void _load(int flag, string inf);
void _initial();
void *_target(void *num);
void _score(string &seq, string &nam, int comp, string flag);

int main(int argc, char* argv[])
{
  _getopt(argc, argv);
  _print("build target database ...");
  _load(1, target);
  _print("load srna sequence ...");
  _load(2, srna);
  _initial();
  otfile.open(output.c_str());
  otfile << "sRNA id\tTarget id\tTarget start\tTarget end\tScore\tFree energy\tDirection" << endl;
  if (rnaSeq.size() < thread)
    thread = rnaSeq.size();
  struct thr tid[thread];
  pthread_mutex_init(&io_mutex, NULL);
  for (int i = 0; i < thread; ++i) {
    _print("target in thread " + _itos(i+1) + " ...");
    tid[i].id = i+1;
    int ret = pthread_create(&tid[i].thread, NULL, _target, &tid[i]);
    if (ret != 0)
      cerr << "pthrea_create error: error_code = " << ret << endl;
  }
  void *status;
  for (int i = 0; i < thread; ++i) {
    int ret = pthread_join(tid[i].thread, &status);
    if (ret != 0)
      cerr << "pthread_join error: error_code = " << ret << endl;
  }
  pthread_mutex_destroy(&io_mutex);
  otfile.close();
  _print("finish all jobs ...");
  return 0;
}

void _getopt(int argc, char* argv[])
{
  char const * shortOpt = "s:t:n:l:m:a:f:N:B:E:P:o:h";
  struct option longOpt[] = {
    {"srna", 1, NULL, 's'},
    {"target", 1, NULL, 't'},
    {"thread", 1, NULL, 'n'},
    {"length", 1, NULL, 'l'},
    {"mirmax", 1, NULL, 'm'}, 
    {"alncut", 1, NULL, 'a'},
    {"fecut", 1, NULL, 'f'},
    {"seednum", 1, NULL, 'N'},
    {"seedbeg", 1, NULL, 'B'},
    {"seedend", 1, NULL, 'E'},
    {"seedpen", 1, NULL, 'P'},
    {"output", 1, NULL, 'o'},
    {"help", 0, NULL, 'h'},
    {NULL, 0, NULL, 0},
  };
  int nextOpt;
  while ((nextOpt = getopt_long(argc, argv, shortOpt, longOpt, NULL)) != -1) {
    switch (nextOpt) {
      case 's':
        srna = optarg;
        break;
      case 't':
        target = optarg;
        break;
      case 'n':
        thread = atoi(optarg);
        break;
      case 'l':
        length = atoi(optarg);
        break;
      case 'm':
        mir_max = atoi(optarg);
        break;
      case 'a':
        score_aln = atof(optarg);
        break;
      case 'f':
        score_fe = atof(optarg);
        break;
      case 'N':
        seed_num = atoi(optarg);
        break;
      case 'B':
        seed_beg = atoi(optarg);
        break;
      case 'E':
        seed_end = atoi(optarg);
        break;
      case 'P':
        seed_pen = atof(optarg);
        break;
      case 'o':
        output = optarg;
        break;
      case 'h':
        help = true;
        break;
    }
  }
  if (output.compare("") == 0)
    output = "sRNATarPredictor_output.xls";
  if (help == true or srna.compare("") == 0 or target.compare("") == 0)
    _usage();
}

void _usage()
{
  printf ("\nqTar: [options]\n");
  printf ("      -s, --srna   *<s>  small rna sequence in fasta format\n");
  printf ("      -t, --target *<s>  target sequence in fasta format\n");
  printf ("      -n, --thread  <i>  number of threads, default 1\n");
  printf ("      -l, --length  <i>  length of k-mer, default 5\n");
  printf ("      -m, --mirmax  <i>  maximum number of mirna target in a target gene, default 5\n");
  printf ("      -a, --alncut  <f>  alignment score threshold, default 18\n");
  printf ("      -f, --fecut   <f>  free energy threshold, default -10\n");
  printf ("      -N, --seednum <i>  minimum number of mirna seed in a candidate target gene, default 2\n");
  printf ("      -B, --seedbeg <i>  start position of mirna seed, default 2\n");
  printf ("      -E, --seedend <i>  end position of mirna seed, default 8\n");
  printf ("      -P, --seedpen <f>  penalty value for mismatch in mirna seed region, default 2\n");
  printf ("      -o, --output  <s>  output, default sRNATarPredictor_output.xls\n");
  printf ("      -h, --help    <b>  print this information\n\n");
  exit(0);
}

void _load(int flag, string inf)
{
  ifstream infile(inf.c_str());
  string line, nam(""), seq;
  int mark = 0;
  while (!infile.eof()) 
  {
    getline(infile, line);
    if (line.find_first_of('>') != string::npos) {
      int end = line.length();
      if (line.find_first_of(' ') != string::npos)
        end = line.find_first_of(' ');
      if (nam.compare("") != 0) {
        transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
        if (flag == 1 && seq.length() > length && tarSeq.find(nam) == tarSeq.end()) {
          tarSeq[nam] = seq;
          for (int i = 0; i < seq.length()-length+1; ++i) {
            tarStr.nam = nam; tarStr.pos = i+1;
            tarPos[seq.substr(i, length)].push_back(tarStr);
          }
        } else if (flag == 2) {
          rnaStr.nam = nam; rnaStr.seq = seq; rnaStr.cmp = _comp(seq);
          rnaSeq[mark % thread == 0 ? thread : mark % thread].push_back(rnaStr);
        }
      }
      nam = line.substr(line.find_first_of('>') + 1, end - line.find_first_of('>') - 1);
      seq = "";
      mark ++;
    } else {
      seq += line;
    }
  }
  transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
  if (flag == 1 && seq.length() > length) {
    tarSeq[nam] = seq;
    for (int i = 0; i < seq.length()-length+1; ++i) {
      tarStr.nam = nam;
      tarStr.pos = i+1;
      tarPos[seq.substr(i, length)].push_back(tarStr);
    }
  } else if (flag == 2) {
    rnaStr.nam = nam; rnaStr.seq = seq; rnaStr.cmp = _comp(seq);
    rnaSeq[mark % thread == 0 ? thread : mark % thread].push_back(rnaStr);
  }
  infile.close();
}

void *_target(void *num)
{
  struct thr *self = (struct thr *) num;
  vector<rnaInf>::iterator itv;
  for (itv = rnaSeq[self->id].begin(); itv != rnaSeq[self->id].end(); ++itv) {
    rnaInf mir = *itv;
    string mir_seq = _revcomp(mir.seq, 0);
    _score(mir_seq, mir.nam, mir.cmp, "forward");
    mir_seq = _revcomp(mir.seq, 1);
    _score(mir_seq, mir.nam, mir.cmp, "reverse");
  }
}

void _score(string &seq, string &nam, int comp, string flag)
{
  map<string, map<int, int> > alnInf; // transcript => [pos => match num]
  for (int i = 0; i < seq.length()-length+1; ++i) {
    string mirSeg = seq.substr(i, length);
    if (tarPos.find(mirSeg) == tarPos.end())
      continue;
    vector<tarInf>::iterator itt;
    for (itt = tarPos[mirSeg].begin(); itt != tarPos[mirSeg].end(); ++itt) {
      tarInf tar = *itt;
      //alnInf[tar.nam + "_" + _itos(tar.pos-i<0 ? 0 : tar.pos-i)] ++;
      alnInf[tar.nam][tar.pos-i<0 ? 0 : tar.pos-i]++;
    }
  }
  for (map<string, map<int, int> >::iterator ita = alnInf.begin(); ita != alnInf.end(); ++ita) {
    vector<PINT> pvec;
    for (map <int, int>::iterator itm = ita->second.begin(); itm != ita->second.end(); ++itm)
      pvec.push_back(make_pair(itm->first, itm->second));
    sort(pvec.begin(), pvec.end(), _pcmp);
    string oline = nam + "\t" + ita->first;
    int oflag = 0;
    for (vector<PINT>::iterator itv = pvec.begin(); itv != pvec.begin() + mir_max; ++itv) {
      if (itv == pvec.end() || itv->second < seed_num)
        break;
      int tend = itv->first + seq.length() > tarSeq[ita->first].length() ? tarSeq[ita->first].length() : itv->first + seq.length();
      string tref = tarSeq[ita->first].substr(itv->first, tend - itv->first + 1);
      // second best match
      int32_t maskLen = strlen(tref.c_str()) / 2;
      maskLen = maskLen < length ? length : maskLen;
      // Declares a default Aligner
      StripedSmithWaterman::Aligner aligner;
      // Declares a default filter
      StripedSmithWaterman::Filter filter;
      // Declares an alignment that stores the result
      StripedSmithWaterman::Alignment alignment;
      aligner.Align(seq.c_str(), tref.c_str(), tref.size(), filter, &alignment, maskLen);
      float saln = alignment.sw_score, sfe = 0.0;
      if (comp == 1)
        sfe += pfe["comp"];
      string pre_base(""), hit_num("");
      int pos_base = 0;
      for (int i = 0; i < alignment.cigar_string.size(); ++i) {
        if (alnCigar.find(alignment.cigar_string[i]) != alnCigar.end()) {
          /* calculate alignment score and free energy */
          int hit = atoi(hit_num.c_str());
          for (int j = 1; j <= hit; ++j) {
            if (pos_base + j >= 2 && pos_base + j <= 8 && alignment.cigar_string[i] != '=')
              saln -= seed_pen;
            char mir_base = seq[pos_base+j-1];
            if (alignment.cigar_string[i] != '=')
              mir_base = 'X';
            if (pre_base.compare("") != 0) {
              pre_base += mir_base;
              if (pfe.find(pre_base) != pfe.end())
                sfe += pfe[pre_base];
            }
            pre_base = mir_base;
            if (hit >= 2 && (alignment.cigar_string[i] == 'I' || alignment.cigar_string[i] == 'D'))
              pre_base = "";
          }
          pos_base += hit;
          hit_num = "";
        } else {
          hit_num += alignment.cigar_string[i];
        }
      }
      if (saln >= score_aln && sfe <= score_fe) {
        oline += "\t" +  _itos(alignment.ref_begin+itv->first+1) +  "\t" + _itos(alignment.ref_end+itv->first+1);
        oline += "\t" + _itos(saln) + "\t" + _itos(sfe) + "\t" + flag;
        oflag ++;
      }
    }
    if (oflag == 0)
      continue;
    pthread_mutex_lock(&io_mutex);
    otfile << oline << endl;
    pthread_mutex_unlock(&io_mutex);
  }
}

/* Freier et al., Improved free-energy parameters for predictions of RNA duplex stability. Proc. Natl. Acad. Sci. USA 83: 9373-9377 (1986)). */
void _initial()
{
  /* Thermodynamic parameters for RNA helix initiation and propagation in 1 M NaCl */
  pfe["init"] = 3.4; pfe["comp"] = 0.4; pfe["ncom"] = 0;
  pfe["AA"] = pfe["TT"] = -0.9; pfe["AT"] = -0.9; pfe["TA"] = -1.1; pfe["CA"] = -1.8;
  pfe["CT"] = -1.7; pfe["GA"] = -2.3; pfe["GT"] = -2.1; pfe["CG"] = -2.0;
  pfe["GC"] = -3.4; pfe["GG"] = pfe["CC"] = -2.9;
  /* Free-energy increments for unpaired terminal nucleotides */
  pfe["AX"] = -0.6; pfe["XA"] = -0.2; pfe["CX"] = -1.2; pfe["XC"] = -0.1;
  pfe["GX"] = -0.6; pfe["XG"] = -0.0; pfe["TX"] = -0.1; pfe["XT"] = -0.2;

  /* alignment cigar */
  alnCigar['='] = alnCigar['X'] = alnCigar['S'] = alnCigar['I'] = alnCigar['D'] = 1;
}

