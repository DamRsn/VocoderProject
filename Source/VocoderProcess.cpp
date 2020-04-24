/*
  ==============================================================================

    VocoderProcess.cpp
    Created: 20 Apr 2020 10:17:42am
    Author:  Damien Ronssin

  ==============================================================================
*/

#include "VocoderProcess.h"

VocoderProcess::VocoderProcess() {    
}

void VocoderProcess::prepare(int wlen_, int hop_, int order_, int orderMax_, std::string windowType) {

    wlen = wlen_;
    hop = hop_;
    startSample = 0;
    window.resize(wlen, 0.0);
    k_iter = 0;
    order = order_;
    orderMax = orderMax_;
    meanAbsE = 34;

    r.resize(orderMax + 1, 1.0);
    a.resize(orderMax + 1, 1.0);
    a_prev.resize(orderMax + 1, 1.0);

    e.resize(wlen, 0.0);
    out.resize(wlen, 0.0);


    if (windowType=="hann")
        dsp::WindowingFunction<float>::fillWindowingTables(&window[0], wlen, dsp::WindowingFunction<float>::hann, false);
    else
    {
        std::cout << "Unknown window type" << std::endl;
        assert(false);
    }
}

VocoderProcess::~VocoderProcess() {}


void VocoderProcess::lpc(const MyBuffer &myBuffer, const int order)
{
    biaisedAutoCorr(myBuffer, r);
    levinsonDurbin(r, a, a_prev);
}


void VocoderProcess::setOrder(int newOrder)
{
    if (newOrder > orderMax)
    {
        std::cerr << "newOrder for LPC is larger than OrderMax" << std::endl;
        assert(false);
    }

    order = newOrder;
}

void VocoderProcess::process(MyBuffer &myBuffer)
{
    while (startSample < myBuffer.getSamplesPerBlock())
    {
        processWindow(myBuffer);
        startSample += hop;
    }

    startSample -= myBuffer.getSamplesPerBlock();
    k_iter+=1;
    k_iter%=200;
}

void VocoderProcess::processWindow(MyBuffer &myBuffer)
{
    int numChannels = myBuffer.getNumChannels();
    float value;
    /*
    for (int channel = 0; channel < numChannels; channel++)
    {
        for (int i = startSample; i < (startSample + wlen); i++)
        {
            value = (myBuffer.getSynthSample(channel, i) + myBuffer.getVoiceSample(channel, i)) * window[i-startSample];
            //value = myBuffer.getVoiceSample(channel, i) * window[i-startSample];
            myBuffer.addOutSample(channel, i, value);
        }
    }
    */
    
    lpc(myBuffer, order);
    filter_FIR(myBuffer, a);
    filter_IIR(myBuffer, a);

}


void VocoderProcess::biaisedAutoCorr(const MyBuffer& myBuffer, std::vector<float>& r)
{
    for (int m = 0; m < order + 1; m++)
    {
        r[m] = 0;

        for (int n = startSample; n < startSample + wlen - m; n++)
        {
            r[m] += myBuffer.getVoiceSample(0, n) * myBuffer.getVoiceSample(0, m + n);
        }
        r[m]/=wlen;
    }
}


void VocoderProcess::levinsonDurbin(const std::vector<float>& r, std::vector<float>& a, std::vector<float>& a_prev)
{
    a[0] = 1.0;
    a[1] = r[1]/r[0];

    double rho_a;
    double r_a;
    double k;

    for (int p = 2;  p < order + 1; p++)
    {
        for (int j = 1; j < p; j++)
        {
            a_prev[j] = a[j];
        }
        rho_a = 0.0;
        r_a = 0.0;

        for (int i = 1; i < p; i++)
        {
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
        a[i]*=-1.;
}


void VocoderProcess::filter_FIR(MyBuffer& myBuffer, const std::vector<float>& a)
{
    meanAbsE = 0.0;
    for (int i = 0; i < wlen; i++) {
        e[i] = a[0] * myBuffer.getVoiceSample(0, startSample + i);

        for (int k = 1; k < order + 1; k++) {
            if (i - k >= 0)
                e[i] += myBuffer.getSynthSample(0, startSample + i - k) * a[k];
            else
                break;
        }
    meanAbsE += abs(e[i]);
    }
    
    meanAbsE/= wlen;
}


void VocoderProcess::filter_IIR(MyBuffer& myBuffer, const std::vector<float>& a) {
    if (k_iter==0)
    {
        std::cout << "meanAbsE:  " << meanAbsE << std::endl;
        std::cout << "a[0] :  " << a[0] << std::endl;
    }
    for (int i = 0; i < wlen; i++) {
        out[i] = myBuffer.getSynthSample(0, startSample + i);
        for (int k = 1; k < order + 1; k++) {
            if (i-k >= 0)
                out[i] -= out[i-k] * a[k];
            else
                break;
        }

        out[i] = out[i]*meanAbsE/(0.2 * a[0]);

        // Todo: something with channels: this code is probably not efficient at all
        // Maybe create a vector out and e for each channel ?
        for (int channel = 0; channel < myBuffer.getNumChannels(); channel++)
            myBuffer.addOutSample(channel, startSample + i, out[i] * window[i]);
    }
}




