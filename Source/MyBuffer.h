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
#include <assert.h>


class MyBuffer
{
public:
    MyBuffer();
    ~MyBuffer();

    void prepare(int samplesPerBlock_, int samplesToKeep_, int latency_, double sampleRate_, int numChannels_);

    void fillInputBuffers(const AudioBuffer<float>& voiceBuffer, const AudioBuffer<float>& synthBuffer);
    void fillOutputBuffer(AudioBuffer<float>& buffer);

    double getVoiceSample(int channel, int idx) const;
    double getSynthSample(int channel, int idx) const;

    void addOutSample(int channel, int idx, float value);
    double getOutSample(int channel, int idx) const;

    int getSamplesPerBlock() const{return samplesPerBlock;}
    int getLatency() const{return latency;}
    int getIdxMax() const {return latency + samplesPerBlock;}
    int getNumChannels() const{return numChannels;}

    void clearOutput(int channel, int numSamples);



private:
    int samplesPerBlock;
    int samplesToKeep;
    double sampleRate;
    int latency;

    int inSize;
    int outSize;
    int numChannels;

    // from which idx to write samplesPerblock samples on mInputVoice and mInputSynth buffers
    int inCounter;

    // from which idx to send samplesPerblock samples in output buffer
    int outCounter;

    // Curr idx corresponding to first sample idx of what will be sent to output at this iteration
    int currCounter;


    AudioBuffer<double> mInputVoice;
    AudioBuffer<double> mInputSynth;
    AudioBuffer<double> mOutput;

};
