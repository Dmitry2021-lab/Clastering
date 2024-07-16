#include <iostream>

#include "LindeBuzoGray.h"
#include "GaussianMixture.h"


int main() {
    TLindeBuzoGray *Clster;
    Clster = new TLindeBuzoGray();

    TGaussianMixture *em;
    em = new TGaussianMixture;
    int nm = 1000;
    int vsize = 11;
    vector<vector<double> > data(nm);

    for (int i = 0; i < nm; i++) {
        data[i].resize(vsize);
        for (int j = 0; j < vsize; j++)
            data[i][j] = rand();
    }
    cout << data[0][2] << endl;
    em->EM_Alg(data, 16);
    Clster->RandomInitPoint(data, 4);
    cout << Clster->LindeBuzoGray_N(data, 16);

}