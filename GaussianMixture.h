#pragma once
#ifndef GaussianMixtureH
#define GaussianMixtureH

#include <numeric>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "stdlib.h"
#include "math.h"

#include "LindeBuzoGray.h"

using namespace std;

const double LSMALL = -0.5e10;
const double EM_THRESH = 1e-4;


class TGaussianMixture
{
private:
	TLindeBuzoGray *lbg;

	double LAdd(double x, double y);
	double PrbState(vector < vector <double> >& SeqCoff, 
		            RModelStage rms, vector <double>& LcPrb, int _time);
	void PstPrbState(vector < vector <double> >& SeqCoff, 
					  RModelStage rms, vector <double>& LcPrb, int _time);
	//RModelStage GaussInit(vector < vector <double> >& SeqCoff, RModelStage InitRms);
public:
	TGaussianMixture();
	~TGaussianMixture();

	RModelStage EM_Alg(vector < vector <double> >& SeqCoff, int NCls);
	double Likelihood(vector < vector <double> >& SeqCoff, RModelStage rms);
};
#endif
