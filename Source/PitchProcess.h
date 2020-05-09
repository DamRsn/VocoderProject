/*
  ==============================================================================

    PitchCorrector.h
    Created: 6 May 2020 10:24:11am
    Author:  Damien Ronssin

  ==============================================================================
*/

#pragma once

#include "../JuceLibraryCode/JuceHeader.h"
#include "MyBuffer.h"
#include "LPC.h"
#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <cassert>


int argExt(MyBuffer& myBuffer, int idxStart, int idxEnd, bool min=true);
class VocoderAudioProcessor;

class PitchProcess
{
public:
    PitchProcess();
    ~PitchProcess();

    void prepare(double fS, double fMin, double fMax, int wlen, int hop, int order, int orderMax, double speed);
    void setAudioProcPtr(VocoderAudioProcessor* audioProcPtr);
    int setLatency();

    void process();

private:

    int yin(MyBuffer& myBuffer);
    void pitchMarks(MyBuffer& myBuffer);
    void getStMarks();
    void psola();

    int wlen;
    int hop;

    double fMin;
    double fMax;
    double fS;
    double delta;

    int tauMax;
    double yinTol;

    int period;
    double pitch;
    double prevPitch;
    double beta;

    // Correction speed in ms
    double speed;

    // Vectors
    // For pitch calculation
    std::vector<double> yinTemp;

    // For pitch marks
    std::vector<int> anMarks;
    std::vector<int> stMarks;
    std::vector<int> anMarksPrev;
    std::vector<int> stMarksPrev;

    // Vectors for LPC
    std::vector<double> a;
    std::vector<double> aPrev;
    std::vector<double> r;
    int order;
    int orderMax;

    // Windows
    std::vector<double> anWindow;
    std::vector<double> stWindow;

    // Pointer to pluginProcessor
    VocoderAudioProcessor* audioProcPtr;

};