/*
  ==============================================================================

    VocoderProcess.cpp
    Created: 20 Apr 2020 10:17:42am
    Author:  Damien Ronssin

  ==============================================================================
*/

#include "VocoderProcess.h"
#include "PluginProcessor.h"
#define PI 3.14159265


VocoderProcess::VocoderProcess(){}

void VocoderProcess::setAudioProcPtr(VocoderAudioProcessor* audioProcPtr)
{
    this->audioProcPtr = audioProcPtr;
}

void VocoderProcess::prepare(int wlen_, int hop_, int order_, int orderMax_, std::string windowType) {
    std::cout.precision(3);
    wlen = wlen_;
    hop = hop_;
    startSample = 0;
    window.resize(wlen, 0.0);
    k_iter = 0;
    orderVoice = order_;
    orderMaxVoice = orderMax_;

    orderMaxSynth = 8;
    orderSynth = 4;

    rVoice.resize(orderMaxVoice + 1, 1.0);
    aVoice.resize(orderMaxVoice + 1, 1.0);
    aPrevVoice.resize(orderMaxVoice + 1, 1.0);
    eVoice.resize(wlen, 0.0);

    rSynth.resize(orderMaxSynth + 1, 1.0);
    aSynth.resize(orderMaxSynth + 1, 1.0);
    aPrevSynth.resize(orderMaxSynth + 1, 1.0);
    eSynth.resize(wlen, 0.0);

    out.resize(wlen, 0.0);

    // gain business
    lambda = 0.0;
    g = 0.0;
    EeSynth = 1.0;
    EeVoice = 0.0;


    // Todo: Change this, put into function and deal with analysis/synthesis windows automatically
    float overlap = float((wlen-hop))/float(wlen);

    float overlapFactor = 1.0;

    if (abs(overlap - 0.75)< pow(10, -10))
    {
        overlapFactor = 1/sqrt(2);
    }

    if (abs(overlap - 0.75) > pow(10, -10) &&
        abs(overlap - 0.5) > pow(10, -10))
    {
        std::cout << "Invalid overlap" << std::endl;
        assert(false);
    }

    if (windowType=="hann")
    {
        dsp::WindowingFunction<float>::fillWindowingTables(&window[0], wlen,
                dsp::WindowingFunction<float>::hann,false);
        for (int i = 0; i < wlen; i++)
        {
            window[i] *= overlapFactor;
        }
    }

    else if (windowType=="sine")
    {
        for (int i = 0; i < wlen; i++)
        {
            window[i] = overlapFactor * sin((i+0.5)*PI/float(wlen));
        }
    }

    else
    {
        std::cout << "Unknown window type" << std::endl;
        assert(false);
    }
}

VocoderProcess::~VocoderProcess() {}


void VocoderProcess::setOrderVoice(int newOrder)
{
    if (newOrder > orderMaxVoice)
    {
        std::cout << "newOrder for LPC is larger than OrderMax" << std::endl;
        assert(false);
    }

    orderVoice = newOrder;
}


void VocoderProcess::setOrderSynth(int newOrder)
{
    if (newOrder > orderMaxSynth)
    {
        std::cout << "newOrder for LPC is larger than OrderMax" << std::endl;
        assert(false);
    }

    orderSynth = newOrder;
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
    k_iter%=300;
}

void VocoderProcess::processWindow(MyBuffer &myBuffer)
{
    int numChannels = myBuffer.getNumChannels();
    /*
    float value;
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

    // LPC on Voice
    lpc(myBuffer, &MyBuffer::getVoiceSample, rVoice, aVoice, aPrevVoice, orderVoice);

    // LPC on Synth
    lpc(myBuffer, &MyBuffer::getSynthSample, rSynth, aSynth, aPrevSynth, orderSynth);

    // Get excitation signal on voice, with energy
    filterFIR(myBuffer, &MyBuffer::getVoiceSample, eVoice, aVoice, orderVoice, EeVoice);

    // Get excitation signal on synth, with energy
    filterFIR(myBuffer, &MyBuffer::getSynthSample, eSynth, aSynth, orderSynth, EeSynth);

    // Final: get final cross-synthesis signal: excitation of synth and enveloppe of voice
    filterIIR(myBuffer, eVoice, aVoice, orderVoice);

}


void VocoderProcess::lpc(MyBuffer& myBuffer, float (MyBuffer::*getSample)(int, int) const, std::vector<float>& r,
        std::vector<float>& a, std::vector<float>& a_prev, const int& order)
{
    biaisedAutoCorr(myBuffer, getSample, r, order);
    levinsonDurbin(r, a, a_prev, order);
}


void VocoderProcess::biaisedAutoCorr(MyBuffer& myBuffer, float (MyBuffer::*getSample)(int, int) const,
        std::vector<float> &r,
        const int& order)
{
    for (int m = 0; m < order + 1; m++)
    {
        r[m] = 0;

        for (int n = 0; n <  wlen - m; n++)
        {
            r[m] += (myBuffer.*getSample)(0, startSample + n) *
                    (myBuffer.*getSample)(0, startSample + m + n);
        }
        r[m]/=float(wlen);
    }
}


void VocoderProcess::levinsonDurbin(const std::vector<float>& r, std::vector<float>& a, std::vector<float>& a_prev,
        const int& order)
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

        float rho_a;
        float r_a;
        float k;

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


void VocoderProcess::filterFIR(MyBuffer& myBuffer, float (MyBuffer::*getSample)(int, int) const, std::vector<float>& e,
        const std::vector<float>& a, const int order, double& E)
{
    E = 0.0;
    for (int i = 0; i < wlen; i++)
    {
        e[i] = a[0] * (myBuffer.*getSample)(0, startSample + i);

        for (int k = 1; k < order + 1; k++)
        {
            if (i - k >= 0)
                e[i] += (myBuffer.*getSample)(0, startSample + i - k) * a[k];
            else
                break;
        }
    E += e[i]*e[i];
    }
}


void VocoderProcess::filterIIR(MyBuffer& myBuffer, const std::vector<float>& e, const std::vector<float>& a,
        const int order)
{
    if (EeSynth > pow(10, -2))
        g = lambda * g + (1-lambda) * sqrt(EeVoice/EeSynth);
    else
        g = 0.0;

    for (int i = 0; i < wlen; i++)
    {
        out[i] =  0.5 * g * eSynth[i];

        // Todo: THIS IS WHERE EXPLOSIONS OFTEN OCCUR: IIR FILTER
        for (int k = 1; k < order + 1; k++) {
            if (i-k >= 0)
                out[i] -= out[i-k] * a[k];
            else
                break;
        }


        // Todo: something with channels: this code is probably not efficient at all
        // Maybe create a vector out and e for each channel ?
        for (int channel = 0; channel < myBuffer.getNumChannels(); channel++)
            myBuffer.addOutSample(channel, startSample + i,
                    (out[i] +
                    0.0 * myBuffer.getSynthSample(0, startSample + i) +
                    0.0 * myBuffer.getVoiceSample(0, startSample + i))
                    * window[i]);
    }
}




