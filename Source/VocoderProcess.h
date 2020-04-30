/*
  ==============================================================================

    VocoderProcess.h
    Created: 20 Apr 2020 10:17:42am
    Author:  Damien Ronssin

  ==============================================================================
*/

#pragma once
#include "MyBuffer.h"
#include "../JuceLibraryCode/JuceHeader.h"
#include <string>
#include <vector>
#include <iostream>
#include <assert.h>
#include <cmath>

class VocoderAudioProcessor;

class VocoderProcess
{
public:
    VocoderProcess();
    ~VocoderProcess();

    void prepare(int wlen, int hop, int order, int orderMax, std::string windowType);

    void process(MyBuffer& myBuffer);

    void setOrderVoice(int newOrder);
    void setOrderSynth(int newOrder);


    void setAudioProcPtr(VocoderAudioProcessor* audioProcPtr);


private:
    int wlen;
    int hop;
    int startSample;

    int orderVoice;
    int orderMaxVoice;
    int orderSynth;
    int orderMaxSynth;

    double meanAbsE;
    int k_iter;

    double g;
    double EeVoice;
    double EeSynth;
    double lambda;

    std::vector<float> window;

    std::vector<float> rVoice;
    std::vector<float> aVoice;
    std::vector<float> aPrevVoice;
    std::vector<float> eVoice;

    std::vector<float> rSynth;
    std::vector<float> aSynth;
    std::vector<float> aPrevSynth;
    std::vector<float> eSynth;


    std::vector<float> out;

    void processWindow(MyBuffer& myBuffer);

    void lpc(MyBuffer& myBuffer, float (MyBuffer::*getSample)(int, int) const, std::vector<float>& r,
            std::vector<float>& a, std::vector<float>& a_prev, const int& order);

    void biaisedAutoCorr(MyBuffer& myBuffer, float (MyBuffer::*getSample)(int, int) const, std::vector<float>& r, const
    int&
    order);

    void levinsonDurbin(const std::vector<float>& r, std::vector<float>& a, std::vector<float>& a_prev,
            const int& order);

    void filterFIR(MyBuffer& myBuffer, float (MyBuffer::*getSample)(int, int) const, std::vector<float>& e,
            const std::vector<float>& a, const int order, double& E);

    void filterIIR(MyBuffer& myBuffer, const std::vector<float>& e, const std::vector<float>& a, const int order);

    VocoderAudioProcessor* audioProcPtr;

};

