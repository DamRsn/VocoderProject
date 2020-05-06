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

double sum(const std::vector<double>& v);
double maxAbs(const std::vector<double>& v);

void shift(std::vector<double>& v);

void VocoderProcess::setAudioProcPtr(VocoderAudioProcessor* audioProcPtr)
{
    this->audioProcPtr = audioProcPtr;
}

void VocoderProcess::prepare(int wlen, int hop, int orderVoice, int orderMax, std::string windowType) {
    std::cout.precision(3);
    this->wlen = wlen;
    this->hop = hop;
    startSample = 0;

    k_iter = 0;
    this->orderVoice = orderVoice;
    this->orderMaxVoice = orderMax;

    this->orderMaxSynth = audioProcPtr->orderMaxSynth;
    this->orderSynth =  audioProcPtr->orderSynth;

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


    EeVoiceArr.resize(10, 0.0);
    EeSynthArr.resize(10, 0.0);

    setWindows(windowType);
}

VocoderProcess::~VocoderProcess() {}

void VocoderProcess::setWindows(std::string windowType)
{
    // Build window vector for analysis and synthesis

    anWindow.resize(wlen, 0.0);
    stWindow.resize(wlen, 0.0);

    double overlap = double((wlen-hop))/double(wlen);

    double overlapFactor = 1.0;

    if (abs(overlap - 0.75) < pow(10, -10))
    {
        overlapFactor = 1.0/sqrt(2);
    }

    if (abs(overlap - 0.75) > pow(10, -10) &&
        abs(overlap - 0.5) > pow(10, -10))
    {
        std::cout << "Invalid overlap" << std::endl;
        assert(false);
    }

    if (windowType=="hann")
    {
        std::cout << "Window type is hann" << std::endl;
        dsp::WindowingFunction<double>::fillWindowingTables(&stWindow[0], wlen,
                                                           dsp::WindowingFunction<double>::hann,false);
        for (int i = 0; i < wlen; i++)
        {
            stWindow[i] *= overlapFactor;
        }

        std::fill(anWindow.begin(), anWindow.end(), 1.0);
    }

    else if (windowType=="sine")
    {
        std::cout << "Window type is sine" << std::endl;
        for (int i = 0; i < wlen; i++)
        {
            anWindow[i] = overlapFactor * sin((i+0.5)*PI/double(wlen));
            stWindow[i] = overlapFactor * sin((i+0.5)*PI/double(wlen));
        }
    }

    else
    {
        std::cerr << "Unknown window type" << std::endl;
        assert(false);
    }
}


void VocoderProcess::setOrderVoice()
{
    if (audioProcPtr->orderVoice > orderMaxVoice)
    {
        std::cout << "New order for LPC Voice is larger than OrderMax" << std::endl;
        assert(false);
    }

    orderVoice = audioProcPtr->orderVoice;
}


void VocoderProcess::setOrderSynth()
{
    if (audioProcPtr->orderSynth > orderMaxSynth)
    {
        std::cout << "New order for LPC Synth is larger than OrderMax" << std::endl;
        assert(false);
    }

    orderSynth = audioProcPtr->orderSynth;
}

void VocoderProcess::process(MyBuffer &myBuffer)
{
    // Get slider values
    gainBeforeIIR = Decibels::decibelsToGain(audioProcPtr->gain);
    gainVoice = Decibels::decibelsToGain(audioProcPtr->gainVoice);
    gainSynth = Decibels::decibelsToGain(audioProcPtr-> gainSynth);
    gainVocoder = Decibels::decibelsToGain(audioProcPtr->gainVocoder);

    // Call process window
    while (startSample < myBuffer.getSamplesPerBlock())
    {
        processWindow(myBuffer);
        startSample += hop;
    }

    // Update startSample
    startSample -= myBuffer.getSamplesPerBlock();

    // Variable used to print something every ... iterations
    k_iter+=1;
    k_iter%=300;
}


void VocoderProcess::processWindow(MyBuffer &myBuffer)
{
    // Update order of LPC
    setOrderSynth();
    setOrderVoice();

    // LPC on Voice
    lpc(myBuffer, &MyBuffer::getVoiceSample, rVoice, aVoice, aPrevVoice, orderVoice);

    // LPC on Synth
    lpc(myBuffer, &MyBuffer::getSynthSample, rSynth, aSynth, aPrevSynth, orderSynth);

    // Get excitation signal on voice, with energy
    filterFIR(myBuffer, &MyBuffer::getVoiceSample, eVoice, aVoice, orderVoice, EeVoice);

    // Get excitation signal on synth, with energy
    filterFIR(myBuffer, &MyBuffer::getSynthSample, eSynth, aSynth, orderSynth, EeSynth);

    // Final: get final cross-synthesis signal: excitation of synth and enveloppe of voice
    filterIIR(myBuffer, aVoice, orderVoice);

}


void VocoderProcess::lpc(MyBuffer& myBuffer, double (MyBuffer::*getSample)(int, int) const, std::vector<double>& r,
        std::vector<double>& a, std::vector<double>& a_prev, const int& order)
{
    biaisedAutoCorr(myBuffer, getSample, r, order);
    levinsonDurbin(r, a, a_prev, order);
}


void VocoderProcess::biaisedAutoCorr(MyBuffer& myBuffer, double (MyBuffer::*getSample)(int, int) const,
        std::vector<double> &r, const int& order)
{
    for (int m = 0; m < order + 1; m++)
    {
        r[m] = 0;

        for (int n = 0; n <  wlen - m; n++)
        {
            r[m] += (myBuffer.*getSample)(0, startSample + n) * anWindow[n] *
                    (myBuffer.*getSample)(0, startSample + m + n) * anWindow[m + n];
        }
        r[m]/=double(wlen);
    }
}


void VocoderProcess::levinsonDurbin(const std::vector<double>& r, std::vector<double>& a, std::vector<double>& a_prev,
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


void VocoderProcess::filterFIR(MyBuffer& myBuffer, double (MyBuffer::*getSample)(int, int) const,
        std::vector<double>& e, const std::vector<double>& a, const int order, double& E)
{
    // Get excitation signal and its energy on current window
    E = 0.0;
    for (int i = 0; i < wlen; i++)
    {
        e[i] = a[0] * (myBuffer.*getSample)(0, startSample + i) * anWindow[i];

        for (int k = 1; k < order + 1; k++)
        {
            if (i - k >= 0)
                e[i] += (myBuffer.*getSample)(0, startSample + i - k) * anWindow[i - k] * a[k];
            else
                break;
        }
    E += e[i]*e[i];
    }
}


void VocoderProcess::filterIIR(MyBuffer& myBuffer, const std::vector<double>& a,
        const int order)
{
    // a is a vector of Voice LPC coefficients
    // order is the order of LPC filter

    // Array of energy of previous windows, use mean to have smooth gain changes
    shift(EeVoiceArr);
    shift(EeSynthArr);
    EeVoiceArr[0] = EeVoice;
    EeSynthArr[0] = EeSynth;


    if (EeSynth > pow(10, -2))
        // For now lambda is 0
        g = lambda * g + (1-lambda) * sqrt(sum(EeVoiceArr)/sum(EeSynthArr));

    else
        g = 0.0;


    for (int i = 0; i < wlen; i++)
    {
        // get first sample and apply gain computed with energy of voice and synth (g)
        // and gain from a slider command (gainBeforeIIR)
        out[i] =  gainBeforeIIR * g * eSynth[i];

        // Todo: THIS IS WHERE EXPLOSIONS OFTEN OCCUR: IIR FILTER
        for (int k = 1; k < order + 1; k++) {
            if (i-k >= 0)
                out[i] -= out[i-k] * a[k];
            else
                break;
        }
    }

    // Add samples to output buffer
    for (int channel = 0; channel < myBuffer.getNumChannels(); channel++)
    {
        for (int i = 0; i < wlen; i++) {
            myBuffer.addOutSample(channel, startSample + i,
                                  (gainVocoder * out[i] +
                                   gainSynth * myBuffer.getSynthSample(0, startSample + i) * anWindow[i] +
                                   gainVoice * myBuffer.getVoiceSample(0, startSample + i) * anWindow[i])
                                  * stWindow[i]);
        }
    }
}


// Some methods on vector
double sum(const std::vector<double>& v)
{
    double s = 0;
    for (auto& n : v)
        s += n;
    return s;
}


double maxAbs(const std::vector<double>& v)
{
    double max = 0.0;
    for (auto& n : v)
    {
        if (abs(n) > max)
            max = abs(n);
    }
    return max;
}


void shift(std::vector<double>& v)
{
    for (int i = v.size()-1; i>0; i--)
    {
        v[i]=v[i-1];
    }
    v[0]=0.0;
}



