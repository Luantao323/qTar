#ifndef GFUN_H_
#define GFUN_H_

#include <vector>
#include <string>
#include <map>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <cstring>
#include <sstream>
#include <algorithm>

using namespace std;

map<char, char> mbase;

string _revcomp(string &str, int flag) 
{
  if (mbase.size() == 0) {
    mbase['A'] = 'T'; mbase['T'] = 'A'; mbase['C'] = 'G'; mbase['G'] = 'C'; mbase['U'] = 'A';
  }
  string revc(str); 
  if (flag == 1)
    reverse(revc.begin(), revc.end());
  for (string::iterator it = revc.begin(); it != revc.end(); ++it)
    if (mbase.find(*it) != mbase.end())
      *it = mbase[*it];
  return revc;
}

void _print(string info)
{
  string str = "[%Y/%m/%d %X] " + info;
  time_t t = time(0);
  char tmp[64];
  strftime (tmp, sizeof(tmp), str.c_str(), localtime(&t));
  puts(tmp);
}

vector<string> _split(string str, char delim)
{
  vector<string> vs;
  string temp;
  for (int i = 0; i < str.length(); ++ i) {
    if (str[i] != delim)
      temp += str[i];
    else {
      if (temp.compare("") != 0)
        vs.push_back(temp);
      temp = "";
    }
  }
  if (temp.compare("") != 0)
    vs.push_back(temp);
  return vs;
}

template <class T>
string _itos(T n)
{
  stringstream ss;
  ss << n;
  return ss.str();
}

int _comp(string &seq) 
{
  int cmp = 0;
  for (int i = 0; i < seq.length() - 2; ++i) {
    string tmp = seq.substr(i, 3);
    if (seq.find(_revcomp(tmp, 1)) != string::npos) {
      cmp = 1;
      break;
    }
  }
  return cmp;
}

#endif
