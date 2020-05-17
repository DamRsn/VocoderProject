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
#include <cassert>
#include <cmath>
#include "add_func.h"
#include "LPC.h"

class VocoderAudioProcessor;

class VocoderProcess
{
public:
    VocoderProcess();
    ~VocoderProcess();

    void prepare(int wlen, int hop, int order, int orderMax, std::string windowType);
    int getLatency(int samplesPerBlock);

    void process(MyBuffer& myBuffer);

    void setAudioProcPtr(VocoderAudioProcessor* audioProcPtr);

private:
    // Window size and hop size, 50% or 75% overlap used
    int wlen;
    int hop;

    // indicator to know where we are in the buffer because wlen, hop and buffer size are different
    int startSample;

    // Order and maximum order for LPC for both Synth and Voice
    int orderVoice;
    int orderMaxVoice;
    int orderSynth;
    int orderMaxSynth;

    // Some stuff to print if needed
    double meanAbsE;
    int k_iter;

    double g;

    // Energy of excitation in current window
    double EeVoice;
    double EeSynth;

    // To update gain factor with an IIR filter of order 1
    double lambda;

    // Different gain to apply
    double gainBeforeIIR;
    double gainVoice;
    double gainVocoder;
    double gainSynth;

    // Analysis and synthesis windows
    std::vector<double> anWindow;
    std::vector<double> stWindow;

    // Autocorr, LPC coeff, Previous LPC coeff (for levinson-durbin), and excitation for voice signal
    std::vector<double> rVoice;
    std::vector<double> aVoice;
    std::vector<double> aPrevVoice;
    std::vector<double> eVoice;

    // Autocorr, LPC coeff, Previous LPC coeff (for levinson-durbin), and excitation for synth signal
    std::vector<double> rSynth;
    std::vector<double> aSynth;
    std::vector<double> aPrevSynth;
    std::vector<double> eSynth;

    // Vocoder output
    std::vector<double> out;

    // Energy of previous and current windows for excitation signal for voice and synth
    std::vector<double> EeVoiceArr;
    std::vector<double> EeSynthArr;


    void setWindows(std::string windowType);

    void setOrderVoice();
    void setOrderSynth();

    void processWindow(MyBuffer& myBuffer);

    void filterFIR(MyBuffer& myBuffer, double (MyBuffer::*getSample)(int, int) const, std::vector<double>& e,
            const std::vector<double>& a, const int order, double& E);

    void filterIIR(MyBuffer& myBuffer, const std::vector<double>& a, const int order);

    // Pointer to pluginProcessor to get the value of the different sliders
    VocoderAudioProcessor* audioProcPtr;

};

