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
#include "Notes.h"
#include <vector>
#include <iostream>
#include <string>
#include <cmath>
#include <cassert>


bool equal(double a, double b);

class VocoderAudioProcessor;

class PitchProcess
{
public:
    PitchProcess();
    ~PitchProcess();

    void prepare(double fS, double fMin, double fMax, int wlen, int hop, int order, int orderMax, double speed);
    void setAudioProcPtr(VocoderAudioProcessor* audioProcPtr);
    int getLatency();

    void process();

private:

    int yin(MyBuffer& myBuffer);
    void pitchMarks(MyBuffer& myBuffer);
    void placeStMarks();
    void psola();

    void fillPsolaWindow(std::vector<double>& psolaWindow, const int& T);
    int argExt(const MyBuffer& myBuffer, int idxStart, int idxEnd, bool min=true);

    int wlen;
    int hop;
    bool valley;

    Notes notes;
    Notes::key key;

    double fMin;
    double fMax;
    double fS;
    double delta;

    int tauMax;
    double yinTol;

    int period;
    int prevPeriod;
    int prevVoicedPeriod;
    int periodNew;
    double pitch;
    double prevPitch;
    double closestFreq;
    double prevClosestFreq;

    double beta;

    // Correction speed in ms
    double speed;
    double counter;

    // Vectors
    // For pitch calculation
    std::vector<double> yinTemp;

    // For pitch marks
    std::vector<int> anMarks;
    std::vector<int> stMarks;
    std::vector<int> prevAnMarks;
    std::vector<int> prevStMarks;

    // Vectors for LPC
    std::vector<double> a;
    std::vector<double> aPrev;
    std::vector<double> r;
    int order;
    int orderMax;

    // Windows
    std::vector<double> anWindow;
    std::vector<double> stWindow;

    std::vector<double> psolaWindow;


    // Pointer to pluginProcessor
    VocoderAudioProcessor* audioProcPtr;

};