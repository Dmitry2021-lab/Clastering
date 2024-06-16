//---------------------------------------------------------------------------
#pragma hdrstop

#include "LindeBuzoGray.h"

//---------------------------------------------------------------------------
#pragma package(smart_init)

TLindeBuzoGray::TLindeBuzoGray()
{
}
//---------------------------------------------------------------------------
TLindeBuzoGray::~TLindeBuzoGray()
{
}
//---------------------------------------------------------------------------
//Take NCls randomly different vectors from array SeqCoff
RModelStage TLindeBuzoGray::RandomInitPoint(vector < vector <double> >& SeqCoff, int NCls)
{
    bool flag;
    int cnt;

    RModelStage Pnt;

    int num = SeqCoff.size();
    int vsize = SeqCoff[0].size();

    vector <int> nPnt(NCls);
    Pnt.Mean.resize(NCls);

    int n = rand() % num;
    cnt = 1;

    nPnt[0] = n;
    Pnt.Mean[0] = SeqCoff[n];

    while (cnt < NCls)
    {
        flag = true;
        n = rand() % num;
        for (int j = 0; j < cnt; j++)
        {
            if (n == nPnt[j])
            {
                flag = false;
                break;
            }
        }
        if (flag == true)
        {
            nPnt[cnt] = n;
            Pnt.Mean[cnt] = SeqCoff[n];
            cnt++;
        }
    }

    return Pnt;
}
//---------------------------------------------------------------------------
double TLindeBuzoGray::EvklDistance(double* Point1, double* Point2, int size)
{
    double sum = 0;
    for (int i = 0; i < size; i++)
        sum += (Point1[i] - Point2[i]) * (Point1[i] - Point2[i]);
    sum = sqrt(sum);
    return sum;
}
//---------------------------------------------------------------------------
RModelStage TLindeBuzoGray::__SetClasterPoint(vector < vector <double> >& SeqCoff, RModelStage Mdl)
{
    double _min, dis;

    int num = SeqCoff.size();
    int vsize = SeqCoff[0].size();
    int cls = -1;

    int limPoints = 0.1*num / Mdl.Cls;

    RModelStage NewMdl;

    vector <int> CntNghPnt(Mdl.Cls, 0);
    vector < vector <double> > Mid(Mdl.Cls, vector <double>(vsize, 0));
    vector < vector <double> > Dsp(Mdl.Cls, vector <double>(vsize, 0));
    vector < vector <double> > lnDsp(Mdl.Cls, vector <double>(vsize, 0));

    for (int i = 0; i < num; i++)
    {
        _min = 1e10;
        for (int j = 0; j < Mdl.Cls; j++)
        {
            dis = EvklDistance(&Mdl.Mean[j][0], &SeqCoff[i][0], vsize);
            if (dis < _min)
            {
                _min = dis;
                cls = j;
            }
        }
        CntNghPnt[cls]++;

        for (int j = 0; j < vsize; j++)
        {
            Mid[cls][j] += SeqCoff[i][j];
            Dsp[cls][j] += SeqCoff[i][j] * SeqCoff[i][j];
        }
    }

    NewMdl.Cls = 0;
    double Max = 0;
    int indexCmp = -1;
    int indexCls = -1;
    int RealNum = 0;

    for (int i = 0; i < Mdl.Cls; i++)
    {
        if (CntNghPnt[i] < limPoints)
            continue;
        RealNum += CntNghPnt[i];
        for (int j = 0; j < vsize; j++)
        {
            Mid[i][j] /= CntNghPnt[i];
            Dsp[i][j] = sqrt(Dsp[i][j] / CntNghPnt[i] - Mid[i][j] * Mid[i][j]);
            lnDsp[i][j] = log(Dsp[i][j]);
            if (Max < Dsp[i][j])
            {
                Max = Dsp[i][j];
                indexCmp = j;
                indexCls = i;
            }
        }
        NewMdl.Mean.push_back(Mid[i]);
        NewMdl.Var.push_back(Dsp[i]);
        NewMdl.lnVar.push_back(lnDsp[i]);
        NewMdl.Cls++;
    }

    for (int i = 0; i < Mdl.Cls; i++)
        NewMdl.lnGaussPrb.push_back(log(CntNghPnt[i]/(double)RealNum));

    NewMdl.indexCmp = indexCmp;
    NewMdl.indexCls = indexCls;
    NewMdl.maxVar = sqrt(Max);

    return NewMdl;
}

//���������� �������������� �������� ������ ��������� � ���������� ������ ������� �����. ��������
vector < vector <int> > TLindeBuzoGray::SetClasterPoint(vector < vector <double> >& SeqCoff, RModelStage Mdl)
{
    double _min = 1e20;
    int index = -1;
    
    double dis;

    int num = SeqCoff.size();
    int vsize = SeqCoff[0].size();

    vector < vector <int> > ClsPoint;
    ClsPoint.resize(Mdl.Cls);

    for (int i = 0; i < num; i++)
    {
        for (int k = 0; k < Mdl.Cls; k++)
        {
            dis = EvklDistance(&Mdl.Mean[k][0], &SeqCoff[i][0], vsize);
            if (dis < _min)
            {
                _min = dis;
                index = k;
            }
        }
        ClsPoint[index].push_back(i);
    }

    //fix ����� ������� � ��������� ����������� �����
    return ClsPoint;
}
//---------------------------------------------------------------------------
RArg TLindeBuzoGray::FirstClaster(vector < vector <double> >& SeqCoff)
{
    RArg DVar;

    int num = SeqCoff.size();
    int vsize = SeqCoff[0].size();

    vector <double> Mean(vsize, 0);
    vector <double> Var(vsize, 0);
    
    double Max = 0;

    for (int j = 0; j < vsize; j++)
    {
        for (int i = 0; i < num; i++)
        {
            Mean[j] += SeqCoff[i][j];
            Var[j] += SeqCoff[i][j] * SeqCoff[i][j];
        }
        Mean[j] /= num;
        Var[j] = Var[j] / num - Mean[j] * Mean[j];
        if (Var[j] > Max)
        {
            Max = Var[j];
            DVar.maxVar = sqrt(Var[j]);
            DVar.indexCmp = j;
            DVar.indexCls = 0;
        }
    }
    DVar.Mean = Mean;
    DVar.Var = Var;

    return DVar;
}
//---------------------------------------------------------------------------
double TLindeBuzoGray::ClastersDistance(RModelStage InitMdl, RModelStage Mdl)
{
    int i, j, k, index;
    double dis, _min, FullDis = 0;
    bool flag;
    int vsize = Mdl.Mean[0].size();
    vector <int> SetIndex;
    SetIndex.push_back(-1);

    for (i = 0; i < Mdl.Cls; i++)
    {
        _min = 1e20;
        for (j = 0; j < InitMdl.Cls; j++)
        {
            flag = false;
            for (k = 0; k < SetIndex.size(); k++)
            {
                if (j == SetIndex[k])
                {
                    flag = true;
                    break;
                }
            }
            if (flag == true)
                continue;
            dis = EvklDistance(&InitMdl.Mean[j][0], &Mdl.Mean[i][0], vsize);
            if (dis < _min)
            {
                _min = dis;
                index = j;
            }
         }
        SetIndex.push_back(index);
        FullDis += _min;
    }
    return FullDis;
}

RArg TLindeBuzoGray::MaxVariances(vector < vector <double> >& SeqCoff, RModelStage Mdl)
{
    RArg DVar;
    double Mean;
    double Var;
    int i, j, k, n;
    int numPntInCls;

    vector < vector <int> > ClsPoint;

    int num = SeqCoff.size();
    int vsize = SeqCoff[0].size();

    DVar.maxVar = 0;
     
    ClsPoint = SetClasterPoint(SeqCoff, Mdl);
    double Max = 0;
    for (k = 0; k < Mdl.Cls; k++)
    {
        numPntInCls = ClsPoint[k].size();
        for (i = 0; i < vsize; i++)
        {
            Mean = 0;
            Var = 0;
            for (j = 0; j < numPntInCls; j++)
            {
                int n = ClsPoint[k][j];
                Mean += SeqCoff[n][i];
                Var += SeqCoff[n][i] * SeqCoff[n][i];
            }
           Mean /= num;
           Var = Var / num - Mean * Mean;
           if (Max < Var)
           {
              Max = Var;
              DVar.maxVar = sqrt(Var);
              DVar.indexCmp = i;
              DVar.indexCls = k;
           }
        }
    }
    return DVar;
}
//---------------------------------------------------------------------------
RModelStage TLindeBuzoGray::K_MeanGeneral(vector <vector <double> >&SeqCoff, RModelStage InitMdl)
{
    double FullDis = 1e10;
    int iterr = 0;

    RModelStage Mdl;
    
    while ((FullDis > KMEAN_THRESH) && (iterr < MAXITERRATION))
    {
        Mdl= __SetClasterPoint(SeqCoff, InitMdl);
        FullDis = ClastersDistance(InitMdl, Mdl);
        InitMdl = Mdl;
        iterr++;
    }
    
    return Mdl;
}

//---------------------------------------------------------------------------
//--------------------------------------------------------------------------- 
//����������� ��������� �����, ���� � ����.
//��������� ��������� �� ������ ��� ��� ���� ��������� ����������
//����� ���� ��������� ���� ���������� ��������� ������������
//�� ��������� ����������� ������������ ���������� �������� ���� �� ���������� ��������������� ����. ���������
//��������� ���������� ����� �������������� ���������� �-mean
RModelStage TLindeBuzoGray::LindeBuzoGray_N(vector < vector <double> >& SeqCoff, int NCls)
{
    double Mult = 0.01;
    
    int num = SeqCoff.size();
    int vsize = SeqCoff[0].size();

    RModelStage Mdl;

    double reliability = 1.0 - 1.0 / sqrt(num / NCls);
    if (reliability < 0.7)
    {
        Mdl.Flag = false;
        cout << "A reliability is not enough!" << endl;
        cout << "You need either more points or less clusters!" << endl;
        return Mdl;
    }

    RArg fCls = FirstClaster(SeqCoff);
    Mdl.Mean.push_back(fCls.Mean);
    if (NCls == 1)
    {
        Mdl.Flag = true;
        Mdl.Cls = NCls;
        return Mdl;
    }


    Mdl.Mean.push_back(fCls.Mean);
    Mdl.Mean[1][fCls.indexCmp] -= fCls.maxVar * Mult;
    Mdl.Mean[0][fCls.indexCmp] += fCls.maxVar * Mult;
    Mdl.Cls = 2;
   

    Mdl = K_MeanGeneral(SeqCoff, Mdl);
    while (Mdl.Cls < NCls)
    {
        fCls = MaxVariances(SeqCoff, Mdl);
        Mdl.Mean.push_back(Mdl.Mean[fCls.indexCls]);
        Mdl.Mean.push_back(Mdl.Mean[fCls.indexCls]);
        int size = Mdl.Mean.size();
        Mdl.Mean[size - 1][fCls.indexCmp] -= fCls.maxVar * Mult;
        Mdl.Mean[size - 2][fCls.indexCmp] += fCls.maxVar * Mult;
        Mdl.Cls += 2;
        Mdl = K_MeanGeneral(SeqCoff, Mdl);
    }
    return Mdl;
}

