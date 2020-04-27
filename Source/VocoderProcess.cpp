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
    this->wlen = wlen_;
    hop = hop_;
    startSample = 0;
    window.resize(wlen, 0.0);
    k_iter = 0;
    order = order_;
    orderMax = orderMax_;
    meanAbsE = 0.0;

    r.resize(orderMax + 1, 1.0);
    a.resize(orderMax + 1, 1.0);
    a_prev.resize(orderMax + 1, 1.0);

    e.resize(wlen, 0.0);
    out.resize(wlen, 0.0);

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


void VocoderProcess::setOrder(int newOrder)
{
    if (newOrder > orderMax)
    {
        std::cout << "newOrder for LPC is larger than OrderMax" << std::endl;
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
    k_iter%=1000000;
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
    
    lpc(myBuffer, order);
    filter_FIR(myBuffer, a, order);
    filter_IIR(myBuffer, a, order);

}


void VocoderProcess::lpc(const MyBuffer &myBuffer, const int order)
{
    biaisedAutoCorr(myBuffer, r, order);
    levinsonDurbin(r, a, a_prev, order);
}


void VocoderProcess::biaisedAutoCorr(const MyBuffer& myBuffer, std::vector<float>& r, const int order)
{
    for (int m = 0; m < order + 1; m++)
    {
        r[m] = 0;

        for (int n = 0; n <  wlen - m; n++)
        {
            r[m] += myBuffer.getVoiceSample(0, startSample + n) *
                    myBuffer.getVoiceSample(0, startSample + m + n);
        }
        r[m]/=float(wlen);
    }
}


void VocoderProcess::levinsonDurbin(const std::vector<float>& r, std::vector<float>& a, std::vector<float>& a_prev,
        const int order)
{
    a[0] = 1.0;
    a[1] = r[1]/r[0];

    float rho_a;
    float r_a;
    float k;

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


void VocoderProcess::filter_FIR(MyBuffer& myBuffer, const std::vector<float>& a, const int order)
{
    meanAbsE = 0.0;
    for (int i = 0; i < wlen; i++) {
        e[i] = a[0] * myBuffer.getVoiceSample(0, startSample + i);

        for (int k = 1; k < order + 1; k++) {
            if (i - k >= 0)
                e[i] += myBuffer.getVoiceSample(0, startSample + i - k) * a[k];
            else
                break;
        }
    meanAbsE += abs(e[i]);
    }
    
    meanAbsE/= wlen;
}


void VocoderProcess::filter_IIR(MyBuffer& myBuffer, const std::vector<float>& a, const int order) {
    int printIdx = 0;

    for (int i = 0; i < wlen; i++)
    {
        out[i] = myBuffer.getSynthSample(0, startSample + i);

        out[i] = out[i]*meanAbsE;


        /*
        if (k_iter== 0)
        {
            std::cout<< out[i] << " ";
            if (i==wlen-1)
                std::cout << "\n \n";
        }
         */

        for (int k = 1; k < order + 1; k++) {
            if (i-k >= 0)
                out[i] -= out[i-k] * a[k];
            else
                break;
        }

        /*
        if(abs(out[i]) > 1. && i > 0)
        {
            out[i] = 0;
            std::cout<< i << "\n";
        }
        */
        /*
        if(abs(out[i])>10)
        {
            std::cout << "out[i]:  " << out[i] << "  meanAbsE:  " << meanAbsE << " \n";

            if (printIdx==0) {
                for (int j = 0; j < order + 1; j++)
                    std::cout << a[j] << "  ";
                std::cout << "\n\n";
                std::cout << "i " << i << "\n";
                std::cout << "x = [";
                for (int j = 0; j < wlen; j++) {
                    std::cout << myBuffer.getSynthSample(0, startSample + j) << ", ";
                }
                std::cout << "] \n";
            }
            printIdx +=1;
        }
        */

        // Todo: something with channels: this code is probably not efficient at all
        // Maybe create a vector out and e for each channel ?
        for (int channel = 0; channel < myBuffer.getNumChannels(); channel++)
            myBuffer.addOutSample(channel, startSample + i,  audioProcPtr->gain * out[i] * window[i]);
    }
}




