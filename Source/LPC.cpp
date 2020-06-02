/*
  ==============================================================================

    LPC.cpp
    Created: 6 May 2020 11:49:54am
    Author:  Damien Ronssin

  ==============================================================================
*/

#include "LPC.h"


void lpc(MyBuffer& myBuffer, double (MyBuffer::*getSample)(int, int) const, std::vector<double>& r,
                         std::vector<double>& a, std::vector<double>& a_prev, const int& order, const int& wlen,
                         const int& startSample, const std::vector<double>& anWindow)
{
    biaisedAutoCorr(myBuffer, getSample, r, order, wlen, startSample, anWindow);
    levinsonDurbin(r, a, a_prev, order);
}


void biaisedAutoCorr(MyBuffer& myBuffer, double (MyBuffer::*getSample)(int, int) const, std::vector<double> &r,
                    const int& order,  const int& wlen, const int& startSample, const std::vector<double>& anWindow)
{
    std::fill(r.begin(), r.end(), 0.0);

    double tmp;
    for (int n=0; n < wlen; n++)
    {
        tmp = (myBuffer.*getSample)(0, startSample + n) * anWindow[n];

        for (int m = 0; m < order + 1; m++)
        {
            if (n >= wlen-m)
                break;
            else
                r[m] +=  tmp * (myBuffer.*getSample)(0, startSample + m + n) * anWindow[m + n];
        }
    }

    for (int m = 0; m < order + 1; m++)
    {
        r[m]/=double(wlen);
    }
}


void levinsonDurbin(const std::vector<double>& r, std::vector<double>& a, std::vector<double>& a_prev, const int& order)
{
    // If first entry of autocorr is 0, then input is full of 0, put a[i] = 0 for all i but 0
    if (abs(r[0])<pow(10, -9))
    {
        std::fill(a.begin(), a.end(), 0.0);
        a[0] = 1.0;
    }

    else {
        a[0] = 1.0;
        a[1] = r[1] / r[0];

        double rho_a;
        double r_a;
        double k;

        for (int p = 2; p < order + 1; p++) {
            for (int j = 1; j < p; j++) {
                a_prev[j] = a[j];
            }
            rho_a = 0.0;
            r_a = 0.0;

            for (int i = 1; i < p; i++) {
                rho_a += r[p - i] * a[i];
                r_a += r[i] * a[i];
            }

            k = (r[p] - rho_a) / (r[0] - r_a);

            for (int i = 1; i < p; i++)
                a[i] = a_prev[i] - k * a_prev[p - i];

            a[p] = k;
        }

        // Multiply all the entries but the first by -1
        for (int i = 1; i < order + 1; i++)
            a[i] *= -1.;
    }
}
