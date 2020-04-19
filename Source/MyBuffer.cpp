/*
  ==============================================================================

    MyBuffer.cpp
    Created: 19 Apr 2020 3:22:24pm
    Author:  Damien Ronssin

  ==============================================================================
*/

#include "MyBuffer.h"

MyBuffer::MyBuffer() {

}

MyBuffer::~MyBuffer() {

}

void MyBuffer::prepare(int samplesPerBlock_, int samplesToKeep_, int latency_, double sampleRate_, int numChannels_)
{
    samplesPerBlock = samplesPerBlock_;
    samplesToKeep = samplesToKeep_;
    sampleRate = sampleRate_;
    latency = latency_;
    numChannels = numChannels_;

    inSize = samplesToKeep + samplesPerBlock + latency;
    outSize = samplesPerBlock + latency;


    mInputVoice.setSize(numChannels, inSize);
    mInputSynth.setSize(numChannels, inSize);

    mOutput.setSize(numChannels, outSize);

    mInputVoice.clear();
    mInputSynth.clear();
    mOutput.clear();

    inCounter = 0;
    outCounter = 0;

    currCounter = (inCounter - latency)%inSize;

}


void MyBuffer::fillInputBuffers(const AudioBuffer<float> &voiceBuffer, const AudioBuffer<float> &synthBuffer)
{
    for(int i = 0; i < numChannels; i++)
    {
        mInputVoice.copyFrom(i, inCounter, voiceBuffer, i, 0, samplesPerBlock);
        mInputSynth.copyFrom(i, inCounter, synthBuffer, i, 0, samplesPerBlock);
    }

    inCounter = (inCounter + samplesPerBlock)%inSize;
    currCounter = (currCounter + samplesPerBlock)%inSize;

}


void MyBuffer::fillOutputBuffer(AudioBuffer<float> &buffer)
{
    buffer.clear();

    for(int i = 0; i < numChannels; i++)
    {
        buffer.copyFrom(i, 0, mOutput, i, outCounter, samplesPerBlock);
    }

    // Clear stuffs that have just been written to buffer
    mOutput.clear(outCounter, samplesPerBlock);

    // update outCounter
    outCounter = (outCounter + samplesPerBlock)%outSize;

}


float MyBuffer::getVoiceSample(int channel, int idx)
{
    if (idx < -samplesToKeep || idx >= (samplesPerBlock + latency))
    {
        std::cerr << "idx out of range in getVoiceSample()." << std::endl;
    }

    return mInputVoice.getSample(channel, (currCounter + idx)%inSize);


}


float MyBuffer::getSynthSample(int channel, int idx)
{
    if (idx < -samplesToKeep || idx >= (samplesPerBlock + latency))
    {
        std::cerr << "idx out of range in getSynthSample()." << std::endl;
    }

    return mInputSynth.getSample(channel, (currCounter + idx)%inSize);


}


void MyBuffer::addOutSample(int channel, int idx, float value)
{
    if (idx < 0 || idx >= (outSize))
    {
        std::cerr << "idx out of range in addOutSample()." << std::endl;
    }

    mOutput.addSample(channel, (outCounter + idx)%outSize, value);
}