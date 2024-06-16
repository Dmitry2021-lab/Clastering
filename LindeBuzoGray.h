//---------------------------------------------------------------------------
#ifndef ClasteringH
#define ClasteringH
//---------------------------------------------------------------------------
#include <numeric>
#include <vector>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "stdlib.h"
#include "math.h"

using namespace std;

const double KMEAN_THRESH = 1e-4;
const int    MAXITERRATION = 200; 

typedef
struct
{
    int indexCmp;
    int indexCls;
    double maxVar;
    vector <double> Mean;
    vector <double> Var;
} RArg;


typedef
struct
{
    bool Flag;
    int Cls;
    int indexCmp;
    int indexCls;
    double maxVar;
    vector <double> lnGaussPrb;
    vector < vector <double> > Mean;
    vector < vector <double> > Var;
    vector < vector <double> > lnVar;
} RModelStage;


class TLindeBuzoGray
{
  private:
      RArg FirstClaster(vector < vector <double> >& SeqCoff);
      RArg MaxVariances(vector < vector <double> >& SeqCoff, RModelStage Mdl);
      vector < vector <int> > SetClasterPoint(vector < vector <double> >& SeqCoff, RModelStage Mdl);
      RModelStage __SetClasterPoint(vector < vector <double> >& SeqCoff, RModelStage Mdl);
      double ClastersDistance(RModelStage InitMdl, RModelStage Mdl);

  public:
      TLindeBuzoGray();
      ~TLindeBuzoGray();
      double EvklDistance(double* Point1, double* Point2, int size);
      RModelStage RandomInitPoint(vector < vector <double> >& SeqCoff, int NCls);
      RModelStage LindeBuzoGray_N(vector < vector <double> >& SeqCoff, int NCls);
      RModelStage K_MeanGeneral(vector <vector <double> >& SeqCoff, RModelStage InitPrm);
};
#endif
