/*
  ==============================================================================

    MyBuffer.cpp
    Created: 19 Apr 2020 3:22:24pm
    Author:  Damien Ronssin

  ==============================================================================
*/

#include "MyBuffer.h"

MyBuffer::MyBuffer() {}

MyBuffer::~MyBuffer(){}

void MyBuffer::prepare(int samplesPerBlock, int samplesToKeep, int latency, double sampleRate, int
numChannelsInVoice, int numChannelsInSynth, int numChannelsOut)
{
    this->samplesPerBlock = samplesPerBlock;
    this->samplesToKeep = samplesToKeep;
    this->sampleRate = sampleRate;
    this->latency = latency;
    this->numChannelsInVoice = numChannelsInVoice;
    this->numChannelsInSynth = numChannelsInSynth;
    this->numChannelsOut = numChannelsOut;


    inSize = samplesToKeep + samplesPerBlock + latency;
    outSize = samplesPerBlock + latency;

    mInputVoice.setSize(numChannelsInVoice, inSize);
    mInputSynth.setSize(numChannelsInSynth, inSize);

    mOutput.setSize(numChannelsOut, outSize);

    // Set everything to 0
    mInputVoice.clear();
    mInputSynth.clear();
    mOutput.clear();

    inCounter = samplesToKeep + latency;
    outCounter = 0;
    currCounter = samplesToKeep;

    voiceReadPtr = mInputVoice.getReadPointer(0);
}


void MyBuffer::fillInputBuffers(const AudioBuffer<float> &voiceBuffer, const AudioBuffer<float> &synthBuffer)
{
    int numSamples = voiceBuffer.getNumSamples();
    for(int channel = 0; channel < numChannelsInVoice; channel++) {
        auto *voiceReadPtr = voiceBuffer.getReadPointer(channel);
        auto *inVoiceWrtPtr = mInputVoice.getWritePointer(channel, 0);

        for (int i = 0; i < numSamples; i++)
            inVoiceWrtPtr[(inCounter + i) % inSize] = voiceReadPtr[i];

    }

    for(int channel = 0; channel < numChannelsInSynth; channel++) {
        auto* synthReadPtr = synthBuffer.getReadPointer(channel);
        auto* inSynthWrtPtr = mInputSynth.getWritePointer(channel, 0);

        if (synthReadPtr!=NULL) {
            for (int i = 0; i < numSamples; i++)
                inSynthWrtPtr[(inCounter + i) % inSize] = synthReadPtr[i];
        }
        else
            for (int i = 0; i < numSamples; i++)
                inSynthWrtPtr[(inCounter + i) % inSize] = 0.0;

    }
}


void MyBuffer::fillOutputBuffer(AudioBuffer<float> &buffer, int nOutputChannels)
{
    buffer.clear();
    for(int channel = 0; channel < nOutputChannels; channel++) {
        auto* channelData = buffer.getWritePointer(channel, 0);
        for (int sample = 0; sample < samplesPerBlock; sample++)
            channelData[sample] = getOutSample(channel, sample);
        
        // Clear stuffs that have just been written to buffer
        clearOutput(channel, samplesPerBlock);
    }


    // update outCounter, inCounter and currCounter
    outCounter = (outCounter + samplesPerBlock)%outSize;
    inCounter = (inCounter + samplesPerBlock)%inSize;
    currCounter = (currCounter + samplesPerBlock)%inSize;
}


double MyBuffer::getVoiceSample(int channel, int idx) const
{
    /*
    if (idx < -samplesToKeep || idx >= (samplesPerBlock + latency)) {
        std::cerr << "idx out of range in getVoiceSample()." << std::endl;
        std::cerr << "idx: " << idx << std::endl;
        assert(false);
    }
     */

    // return mInputVoice.getSample(channel, (currCounter + idx + inSize)%inSize);
    return voiceReadPtr[(currCounter + idx + inSize)%inSize];
}


double MyBuffer::getSynthSample(int channel, int idx) const
{
    /*
    if (idx < -samplesToKeep || idx >= (samplesPerBlock + latency)) {
        std::cerr << "idx out of range in getSynthSample()." << std::endl;
        assert(false);
    }
     */

    return mInputSynth.getSample(channel, (currCounter + idx + inSize)%inSize);
}


void MyBuffer::addOutSample(int channel, int idx, double value)
{
     /*
    if (idx < 0 || idx >= (outSize)) {
        std::cerr << "idx out of range in addOutSample()." << std::endl;
        assert(false);
    }
    */
    mOutput.addSample(channel, (outCounter + idx)%outSize, value);
}


double MyBuffer::getOutSample(int channel, int idx) const
{
    /*
    if (idx < 0 || idx >= (outSize)) {
        std::cerr << "idx out of range in getOutSample()." << std::endl;
        assert(false);
    }
     */
    
    return mOutput.getSample(channel, (outCounter + idx)%outSize);;
}


void MyBuffer::clearOutput(int channel, int numSamples)
{
    if (outCounter + numSamples < outSize)
        mOutput.clear(channel, outCounter, numSamples);

    else
    {
        mOutput.clear(channel, outCounter, outSize-outCounter);
        mOutput.clear(channel, 0, numSamples - outSize + outCounter);
    }
}
