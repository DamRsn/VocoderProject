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
template <typename T>
void operator-=(std::vector<T>& v1, const T& a);
template <typename T>
void operator+=(std::vector<T>& v1, const T& a);
template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v);


class VocoderAudioProcessor;

class PitchProcess
{
public:
    PitchProcess();
    ~PitchProcess();

    void prepare(double fS, double fMin, double fMax, int frameLen, int hop, int order, int orderMax, double speed, int
    samplesPerBlock, Notes::key key);
    void setAudioProcPtr(VocoderAudioProcessor* audioProcPtr);
    int getLatency(int samplesPerBlock);
    int getSampleToKepp();

    //void process(MyBuffer& myBuffer);
    void process(MyBuffer& myBuffer);

private:

    void yin(MyBuffer& myBuffer);
    void pitchMarks(MyBuffer& myBuffer);
    void placeStMarks();
    void psola(MyBuffer& myBuffer);

    /*
    void processStartFrame(MyBuffer& myBuffer);
    void processContFrame(MyBuffer& myBuffer);
    */
    void fillOutputBuffer(MyBuffer& myBuffer);

    void processChunkCont(MyBuffer& myBuffer);
    void processChunkStart(MyBuffer& myBuffer);


    // Utility functions
    void fillPsolaWindow(std::vector<double>& psolaWindow, const int& T);
    int argExt(const MyBuffer& myBuffer, int idxStart, int idxEnd, bool min=true);
    int getClosestAnMarkIdx(const std::vector<int>& anMarks, const std::vector<int>& prevAnMarks, const int& stMark,
            int periodPsola);
    void buildWindows(std::vector<double>& anWindow, std::vector<double>& stWindow);

    // Interpolate
    void interp(std::vector<double>& x, const std::vector<double>& y, std::vector<double>& outWindow,
            const std::vector<double>& psolaWindow, const int& startIdx, const int& stopIdx);

    int frameLen;
    int hop;
    double overlap;
    int startSample;
    int bufferIdxMax;

    int samplesPerBlock;

    int nChunk;
    int chunkSize;
    int chunksPerFrame;


    Notes notes;
    Notes::key key;

    bool valley;
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
    double prevVoicedPitch;
    double closestFreq;
    double prevClosestFreq;

    int stMarkIdx;
    int nAnMarksOv;
    int nStMarksOv;
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
    std::vector<double> e;
    int order;
    int orderMax;

    // Windows
    std::vector<double> anWindow;
    std::vector<double> stWindow;

    // vector with hann window for psola
    std::vector<double> psolaWindow;

    // to store 2*period + 1 samples centered on a pitch mark
    std::vector<double> periodSamples;

    std::vector<double> xInterp;

    std::vector<double> outFrame;




    // Pointer to pluginProcessor
    VocoderAudioProcessor* audioProcPtr;

};