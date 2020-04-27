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

    void setOrder(int newOrder);

    void setAudioProcPtr(VocoderAudioProcessor* audioProcPtr);


private:
    int wlen;
    int hop;
    int startSample;

    int order;
    int orderMax;
    float meanAbsE;
    int k_iter;
    std::vector<float> window;

    std::vector<float> r;
    std::vector<float> a;
    std::vector<float> a_prev;

    std::vector<float> e;
    std::vector<float> out;

    void processWindow(MyBuffer& myBuffer);

    void biaisedAutoCorr(const MyBuffer& myBuffer, std::vector<float>& r, const int order);
    void levinsonDurbin(const std::vector<float>& r, std::vector<float>& a, std::vector<float>& a_prev,
            const int order);
    void lpc(const MyBuffer& myBuffer, const int order);

    void filter_FIR(MyBuffer& myBuffer, const std::vector<float>& a, const int order);
    void filter_IIR(MyBuffer& myBuffer, const std::vector<float>& a, const int order);

    VocoderAudioProcessor* audioProcPtr;

};

