/*
  ==============================================================================

    MyBuffer.h
    Created: 19 Apr 2020 3:22:24pm
    Author:  Damien Ronssin

  ==============================================================================
*/

#pragma once

#include "../JuceLibraryCode/JuceHeader.h"
#include <iostream>
#include <cmath>
#include <cassert>


class MyBuffer
{
public:
    MyBuffer();
    ~MyBuffer();

    void prepare(int samplesPerBlock, int samplesToKeep, int latency, double sampleRate, int numChannelsInVoice, int
    numChannelsInSynth, int numChannelsOut);

    void fillInputBuffers(const AudioBuffer<float>& voiceBuffer, const AudioBuffer<float>& synthBuffer);
    void fillOutputBuffer(AudioBuffer<float>& buffer, int nOutputChannels);

    double getVoiceSample(int channel, int idx) const;
    double getSynthSample(int channel, int idx) const;

    void addOutSample(int channel, int idx, double value);
    double getOutSample(int channel, int idx) const;

    int getSamplesPerBlock() const {return samplesPerBlock;}
    int getLatency() const {return latency;}
    int getIdxMax() const {return latency + samplesPerBlock;}
    int getSamplesToKeep() const {return samplesToKeep;}
    int getNumOutChannels() const {return numChannelsOut;}
    int getInSize() const {return inSize;}
    int getCurrCounter() const {return currCounter;}

    const double* getVoiceReadPtr() const{return mInputVoice.getReadPointer(0);}
    const double* getSynthReadPtr(int channel) const{return mInputSynth.getReadPointer(channel);}

    double getRMSLevelVoice(int startSample, int numSamples);
    double getRMSLevelSynth(int startSample, int numSamples);

    double getRMSLevelVoiceFull();
    double getRMSLevelSynthFull();

    void addDryVoice(double gain);
    void addSynth(double gain);

    void clearOutput(int channel, int numSamples);

private:
    int samplesPerBlock;
    int samplesToKeep;
    double sampleRate;
    int latency;

    int inSize;
    int outSize;

    int numChannelsInVoice;
    int numChannelsInSynth;
    int numChannelsOut;

    // from which idx to write samplesPerblock samples on mInputVoice and mInputSynth buffers
    int inCounter;

    // from which idx to send samplesPerblock samples in output buffer
    int outCounter;

    // Curr idx corresponding to first sample idx of what will be sent to output at this iteration
    int currCounter;

    // Audio Buffer, implemented circularly in this class
    AudioBuffer<double> mInputVoice;
    AudioBuffer<double> mInputSynth;
    AudioBuffer<double> mOutput;

    const double* voiceReadPtr;

};
