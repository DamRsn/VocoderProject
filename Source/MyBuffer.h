/*
  ==============================================================================

    MyBuffer.h
    Created: 19 Apr 2020 3:22:24pm
    Author:  Damien Ronssin

  ==============================================================================
*/

#pragma once

#include "../JuceLibraryCode/JuceHeader.h"
#include "iostream"


class MyBuffer
{
public:
    MyBuffer();
    ~MyBuffer();

    void prepare(int samplesPerBlock_, int samplesToKeep_, int latency_, double sampleRate_, int numChannels_);

    void fillInputBuffers(const AudioBuffer<float>& voiceBuffer, const AudioBuffer<float>& synthBuffer);
    void fillOutputBuffer(AudioBuffer<float>& buffer);

    float getVoiceSample(int channel, int idx);
    float getSynthSample(int channel, int idx);

    void addOutSample(int channel, int idx, float value);


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


    AudioBuffer<float> mInputVoice;
    AudioBuffer<float> mInputSynth;
    AudioBuffer<float> mOutput;




};