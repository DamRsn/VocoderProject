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

    void prepare(double fS, double fMin, double fMax, int frameLen, int hop, int samplesPerBlock, double
    silenceThresholdDb);
    void prepare2 (MyBuffer& myBuffer);
    void setAudioProcPtr(VocoderAudioProcessor* audioProcPtr);
    int getLatency(int samplesPerBlock);
    int getSampleToKeep();

    void process(MyBuffer& myBuffer);
    void silence();

private:

    void yin(const MyBuffer& myBuffer);
    void computeYinTemp(const MyBuffer& myBuffer, std::vector<double>& yinTemp);
    void pitchMarks(const MyBuffer& myBuffer);
    void placeStMarks();
    void psola(MyBuffer& myBuffer);

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
            const int& startIdx, const int& stopIdx);

    void filterFIR(MyBuffer& myBuffer, int startIdxBuf, int samplesToFilter, int startIdxE);
    void filterIIR();


    int frameLen;
    int hop;
    double overlap;
    int startSample;
    int bufferIdxMax;

    double silenceThresholdDb;

    int samplesPerBlock;

    int nChunk;
    int chunkSize;
    int chunksPerFrame;

    // Value from myBuffer
    int samplesToKeep;

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
    std::vector<double> eFrame;
    int order;
    int orderMax;

    // Windows
    std::vector<double> anWindow;
    std::vector<double> stWindow;

    // vector with hann window for psola
    std::vector<double> psolaWindow;

    // to store 2 * period + 1 samples centered on a pitch mark
    std::vector<double> periodSamples;

    std::vector<double> xInterp;

    std::vector<double> outEFrame;
    std::vector<double> yFrame;

    // Pointer to pluginProcessor
    VocoderAudioProcessor* audioProcPtr;

};