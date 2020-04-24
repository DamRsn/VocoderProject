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

MyBuffer::~MyBuffer()
{
    mOutput.~AudioBuffer();
    mInputSynth.~AudioBuffer();
    mInputVoice.~AudioBuffer();
}

void MyBuffer::prepare(int samplesPerBlock_, int samplesToKeep_, int latency_, double sampleRate_, int numChannels_)
{
    samplesPerBlock = samplesPerBlock_;
    samplesToKeep = samplesToKeep_;
    sampleRate = sampleRate_;
    latency = latency_;
    numChannels = numChannels_;

    if (samplesToKeep != 0)
    {
        std::cerr<< "ERROR: sampleToKeep is not zero, assert call" << std::endl;
        assert(false);
    }
    
    
    inSize = samplesToKeep + samplesPerBlock + latency;
    outSize = samplesPerBlock + latency;


    mInputVoice.setSize(numChannels, inSize);
    mInputSynth.setSize(numChannels, inSize);

    mOutput.setSize(numChannels, outSize);

    // Set everything to 0
    mInputVoice.clear();
    mInputSynth.clear();
    mOutput.clear();

    inCounter = samplesToKeep + latency;
    outCounter = 0;
    currCounter = samplesToKeep;
}


void MyBuffer::fillInputBuffers(const AudioBuffer<float> &voiceBuffer, const AudioBuffer<float> &synthBuffer)
{
    for(int channel = 0; channel < numChannels; channel++)
    {
        if (inCounter + samplesPerBlock <= inSize)
        {
            mInputVoice.copyFrom(channel, inCounter, voiceBuffer, channel, 0, samplesPerBlock);
            mInputSynth.copyFrom(channel, inCounter, synthBuffer, channel, 0, samplesPerBlock);
        }
        
        else
        {
            mInputVoice.copyFrom(channel, inCounter, voiceBuffer, channel, 0, inSize - inCounter);
            mInputSynth.copyFrom(channel, inCounter, synthBuffer, channel, 0, inSize - inCounter);
            
            int n_to_add = samplesPerBlock + inCounter - inSize;
            mInputVoice.copyFrom(channel, 0, voiceBuffer, channel, inSize - inCounter, n_to_add);
            mInputSynth.copyFrom(channel, 0, synthBuffer, channel, inSize - inCounter, n_to_add);
        }
    }
}


void MyBuffer::fillOutputBuffer(AudioBuffer<float> &buffer)
{
    buffer.clear();

    for(int channel = 0; channel < numChannels; channel++)
    {
        auto* channelData = buffer.getWritePointer(channel, 0);
        for (int sample = 0; sample < samplesPerBlock; sample++)
        {
            channelData[sample] = getOutSample(channel, sample);
        }

        // Clear stuffs that have just been written to buffer
        clearOutput(channel, samplesPerBlock);
    }

    // update outCounter, inCounter and currCounter
    outCounter = (outCounter + samplesPerBlock)%outSize;
    inCounter = (inCounter + samplesPerBlock)%inSize;
    currCounter = (currCounter + samplesPerBlock)%inSize;
}


float MyBuffer::getVoiceSample(int channel, int idx) const
{
    if (idx < -samplesToKeep || idx >= (samplesPerBlock + latency))
    {
        std::cerr << "idx out of range in getVoiceSample()." << std::endl;
        std::cerr << "idx: " << idx << std::endl;
        assert(false);
    }

    return mInputVoice.getSample(channel, (currCounter + idx)%inSize);
}


float MyBuffer::getSynthSample(int channel, int idx) const
{
    if (idx < -samplesToKeep || idx >= (samplesPerBlock + latency))
    {
        std::cerr << "idx out of range in getSynthSample()." << std::endl;
        assert(false);
    }

    return mInputSynth.getSample(channel, (currCounter + idx)%inSize);
}


void MyBuffer::addOutSample(int channel, int idx, float value)
{
    if (idx < 0 || idx >= (outSize))
    {
        std::cerr << "idx out of range in addOutSample()." << std::endl;
        assert(false);
    }

    mOutput.addSample(channel, (outCounter + idx)%outSize, value);
}


float MyBuffer::getOutSample(int channel, int idx) const
{
    if (idx < 0 || idx >= (outSize))
    {
        std::cerr << "idx out of range in getOutSample()." << std::endl;
        assert(false);
    }

    return mOutput.getSample(channel, (outCounter + idx)%outSize);
}


void MyBuffer::clearOutput(int channel, int numSamples)
{
    auto* wrtPtr = mOutput.getWritePointer(channel, 0);

    for(int i = 0; i < numSamples; i++)
    {
        wrtPtr[(outCounter + i)%outSize] = 0.0;
    }
}
