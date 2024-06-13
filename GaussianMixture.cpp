#pragma hdrstop

#include "GaussianMixture.h"

#pragma package(smart_init)
//---------------------------------------------------------------------------

TGaussianMixture::TGaussianMixture()
{
    lbg = new TLindeBuzoGray();
}
TGaussianMixture::~TGaussianMixture()
{
    delete lbg;
}
//---------------------------------------------------------------------------
double TGaussianMixture::LAdd(double x, double y)
{

    float minLogExp = -30.0;
    float temp, diff, z;

    if (x < y)
    {
        temp = x;
        x = y;
        y = temp;
    }
    diff = y - x;
    if (diff < minLogExp)
    {
        if (x < LSMALL) return -1e10;
        else return x;
    }
    else
    {
        z = exp(diff);
        return x + log(1.0 + z);
    }
}
/*
RModelStage TGaussianMixture::GaussInit(vector < vector <double> >& SeqCoff, RModelStage InitRms)
{
    int i, j, k, ind;
    double dis, _min, par;
    bool flag;

    int num = SeqCoff.size();
    int vsize = SeqCoff[0].size();

    int limPoints = 0.1 * num / InitRms.Cls;

    RModelStage rms;

    vector < vector <double> > SumCoff(InitRms.Cls, vector <double>(vsize, 0));
    vector < vector <double> > DspSumCoff(InitRms.Cls, vector <double>(vsize, 0));
    vector <int> CntNghPnt(InitRms.Cls, 0);

    for (i = 0; i < num; i++)
    {
        _min = 1e10;
        ind = -1;
        for (j = 0; j < InitRms.Cls; j++)
        {
            dis = lbg->EvklDistance(&InitRms.Mean[j][0], &SeqCoff[i][0], vsize);
            if (dis < _min)
            {
                _min = dis;
                ind = j;
            }
        }

        for (k = 0; k < vsize; k++)
        {
            SumCoff[ind][k] += SeqCoff[i][k];
            DspSumCoff[ind][k] += (SeqCoff[i][k] * SeqCoff[i][k]);
        }
        CntNghPnt[ind]++;
    }

    rms.Cls = 0;
    for (j = 0; j < InitRms.Cls; j++)
    {
        if (CntNghPnt[j] > limPoints)
        {
            flag = true;
            rms.lnAPrioriGaussPrb[rms.CntCls] = log((double)CntNghPnt[j] / (double)num);
            for (k = 0; k < vsize; k++)
            {
                rms.Mean[rms.Cls][k] = SumCoff[j][k] / (double)CntNghPnt[j];
                par = DspSumCoff[j][k] / (double)CntNghPnt[j] - rms.Mid[rms.Cls][k] * rms.Mid[rms.Cls][k];
                if (rms.Mid[rms.Cls][k] != 0 && par <= 0)
                {
                    flag = false;
                    par = 1;
                }
                else if (rms.Mid[rms.Cls][k] == 0)
                {

                    rms.Var[rms.Cls][k] = 0;
                    rms.lnVar[rms.Cls][k] = LSMALL;
                }
                else
                {
                    rms.Var[rms.Cls][k] = sqrt(par);
                    rms.lnVar[rms.Cls][k] = log(rms.Var[rms.Cls][k]);
                }
            }
            if (flag == true) rms.Cls++;
        }
        else
    }

    return rms
}
*/
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
double TGaussianMixture::PrbState(vector < vector <double> >& SeqCoff, RModelStage rms, 
                                   vector <double> &LcPrb, int _time)
{
    double par, Prb;
    double sum = LSMALL;

    int vsize = SeqCoff[0].size();

    for (int j = 0; j < rms.Cls; j++)
    {
        Prb = 0;
        for (int k = 0; k < vsize; k++)
        {
            par = (rms.Mean[j][k] - SeqCoff[_time][k]) / rms.Var[j][k];
            Prb -= rms.lnVar[j][k] + par * par / 2.0;
        }
        LcPrb[j] = rms.lnGaussPrb[j] + Prb;
        sum = LAdd(sum, LcPrb[j]);
    }

    return sum;
}
//---------------------------------------------------------------------------
void TGaussianMixture::PstPrbState(vector < vector <double> >& SeqCoff, RModelStage rms,
                                    vector <double>& LcPrb, int _time)
{
    int i, j, k;
    double sum, Prb;
    long double par;
    int vsize = SeqCoff[0].size();

    sum = PrbState(SeqCoff, rms, LcPrb, _time);

    for (j = 0; j < rms.Cls; j++)
        LcPrb[j] = exp(LcPrb[j] - sum);

}
//---------------------------------------------------------------------------
double TGaussianMixture::Likelihood(vector < vector <double> >& SeqCoff, RModelStage rms)
{
    vector <double> LcPrb(rms.Cls, 0);
    double sum, sum1 = LSMALL;

    int num = SeqCoff.size();
    int vsize = SeqCoff[0].size();

    for (int _time = 0; _time < num; _time++)
    {
        sum = PrbState(SeqCoff, rms, LcPrb, _time);
        sum1 = LAdd(sum1, sum);
    }
    return sum1;
}
//---------------------------------------------------------------------------
RModelStage TGaussianMixture::EM_Alg(vector < vector <double> >& SetVectors, int NCls)
{
    bool Flag;
    int n, i, j;
    double StopIt;
    double sum1, sum2, par_d, par_m, par_ld;
    double lh, Oldlh;
    int rcs = 0;
    
    RModelStage rms, NewRms;

    
    int num = SetVectors.size();
    int vsize = SetVectors[0].size();
 
    NewRms.Cls = rms.Cls;
    
    StopIt = 100;
    lh = -1000;

    rms = lbg->LindeBuzoGray_N(SetVectors, NCls);
    
    vector < vector <double> > LockPrb(num, vector <double>(rms.Cls, 0));

    while ((StopIt > EM_THRESH) && (rcs < MAXITERRATION))
    {
        
        vector <double> sum(rms.Cls);
        for (n = 0; n < num; n++)
        {
            PstPrbState(SetVectors, rms, LockPrb[n], n);
            for (j = 0; j < rms.Cls; j++)
                sum[j] += LockPrb[n][j];
        }

        for (j = 0; j < rms.Cls; j++)
            NewRms.lnGaussPrb.push_back(log(sum[j] / (double)num));

        NewRms.Cls = 0;
        for (j = 0; j < rms.Cls; j++)
        {
            Flag = true;
            vector <double> current_m;
            vector <double> current_d;
            vector <double> current_ld;
            for (i = 0; i < vsize; i++)
            {
                sum1 = 0;
                sum2 = 0;
                for (int k = 0; k < num; k++)
                {
                    sum1 += SetVectors[k][i] * LockPrb[k][j];
                    sum2 += SetVectors[k][i] * SetVectors[k][i] * LockPrb[k][j];
                }
                par_m = sum1 / sum[j];
                par_d = sqrt(sum2 / sum[j] - par_m * par_m);
                par_ld = log(par_d);
                if (par_d > 1e-7)
                {
                    current_m.push_back(par_m);
                    current_d.push_back(par_d);
                    current_ld.push_back(par_ld);
                }
                else
                {
                    Flag = false;
                    break;
                }

            }
            if (Flag == true)
            {
                NewRms.Mean.push_back(current_m);
                NewRms.Var.push_back(current_d);
                NewRms.lnVar.push_back(current_ld);
                NewRms.Cls++;
            }
        }

        if (NewRms.Cls == rms.Cls)
        {
            Oldlh = lh;
            lh = Likelihood(SetVectors, rms);
            StopIt = (lh - Oldlh);
        }
        else
        {
            StopIt = 100;
        }
        rms = NewRms;
        NewRms.Cls = 0;
        NewRms.lnGaussPrb.clear();
        NewRms.lnVar.clear();
        NewRms.Var.clear();
        NewRms.Mean.clear();
        rcs++;
    }
    return rms;
}
//---------------------------------------------------------------------------
