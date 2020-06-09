/*
  ==============================================================================

    LPC.cpp
    Created: 6 May 2020 11:49:54am
    Author:  Damien Ronssin

  ==============================================================================
*/

#include "LPC.h"

/**
 * Compute LPC and update the LPC coefficients in vector a
 * @param myBuffer : reference to my buffer
 * @param inputBuffer : pointer to inputBuffer (either voice or synth). Use this for vectorization
 * @param r : vector for autocorrelation
 * @param a : vector for filter coeff
 * @param a_prev : vector for filter coeff calculation
 * @param order : order of LPC filter
 * @param wlen : length of a frame
 * @param startSample : needed to access values in myBuffer at the right place
 * @param anWindow : analysis window (vector of size wlen)
 */
void lpc(MyBuffer& myBuffer, const double* inputBuffer/*double (MyBuffer::*getSample)(int, int) const*/,
        std::vector<double>& r,
                         std::vector<double>& a, std::vector<double>& a_prev, const int& order, const int& wlen,
                         const int& startSample, const std::vector<double>& anWindow)
{
    biaisedAutoCorr(myBuffer, inputBuffer, r, order, wlen, startSample, anWindow);
    levinsonDurbin(r, a, a_prev, order);
}

/**
 * Compute the biased autocorrelation needed for LPC calculation in Levinson-Durbin recursion
 * @param myBuffer
 * @param inputBuffer
 * @param r
 * @param order
 * @param wlen
 * @param startSample
 * @param anWindow
 */
void biaisedAutoCorr(MyBuffer& myBuffer, const double* inputBuffer/*double (MyBuffer::*getSample)(int, int) const*/,
        std::vector<double> &r,
                    const int& order,  const int& wlen, const int& startSample, const std::vector<double>& anWindow)
{
    std::fill(r.begin(), r.end(), 0.0);

    // idx for myBuffer
    int currCounter = myBuffer.getCurrCounter();
    int inSize = myBuffer.getInSize();

    int startIdx;
    int stopIdx;

    double tmp;

    for (int n=0; n < wlen; n++)
    {
        tmp = inputBuffer[(currCounter + startSample + n + inSize)%inSize] * anWindow[n];

        startIdx = (currCounter + startSample + n + inSize)%inSize;
        stopIdx = (currCounter + startSample + n + order + inSize)%inSize;

        if (startIdx < stopIdx) {
            #pragma clang loop vectorize(enable)
            for (int m = 0; m < order + 1; m++) {
                if (n < wlen - m)
                    r[m] += tmp * inputBuffer[startIdx + m] * anWindow[m + n];
            }
        }
        else
        {
            int k = inSize - startIdx;

            #pragma clang loop vectorize(enable)
            for (int m = 0; m < k; m++) {
                if (n < wlen - m)
                    r[m] += tmp * inputBuffer[startIdx + m] * anWindow[m + n];
            }

            #pragma clang loop vectorize(enable)
            for (int m = k; m < order + 1; m++)
            {
                if (n < wlen - m)
                    r[m] += tmp * inputBuffer[m - k] * anWindow[m + n];
            }
        }
    }

    #pragma clang loop vectorize(enable)
    for (int m = 0; m < order + 1; m++)
    {
        r[m]/=double(wlen);
    }
}


/**
 * Levinson-Durbin Recursion. Uptdates the vector a (filter coefficients) given the autocorrelation r
 * @param r
 * @param a
 * @param a_prev
 * @param order
 */
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
