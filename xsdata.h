#ifndef XSDATA_H
#define XSDATA_H

#include <string>
#include <vector>

using namespace std;

static vector<double> xs_E;
static vector<double> xs_v;
static vector<double> xs_sig;

int readxs(const string &filename, vector<double> &xs_E, vector<double> &xs_sig);
void gridEtoV(vector<double> &xs_E, vector<double> &xs_v);

#endif
