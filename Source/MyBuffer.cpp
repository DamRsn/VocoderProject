/*
  ==============================================================================

    MyBuffer.cpp
    Created: 19 Apr 2020 3:22:24pm
    Author:  Damien Ronssin

  ==============================================================================
*/

#include "MyBuffer.h"

/**
 * Class with circular buffer of double precision numbers for voice (mono), synthesizer (stereo) and output(stereo)
 * signals
 */

MyBuffer::MyBuffer() {}

MyBuffer::~MyBuffer(){}


/**
 * Function to call in prepare to play to allocate the memory of the different buffers and initialised the
 * different variables
 * @param samplesPerBlock : buffer size
 * @param samplesToKeep : number of past samples to keep (to be accessed with negative number in getSample methods)
 * @param latency : number of samples of latency
 * @param sampleRate : in Hz
 * @param numChannelsInVoice : number of channels for voice (1)
 * @param numChannelsInSynth : number of channels for synthesizer (2)
 * @param numChannelsOut : number of channel for output buffer (2)
 */
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

    // Size of voice and synthesizer buffers (in number of samples)
    inSize = samplesToKeep + samplesPerBlock + latency;
    // Size of output buffer in number of samples
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


/**
 * Replace last samplesPerBlock samples by samplesPerBlock samples from current buffer to mInputVoice and mInputSynth
 * cast all the sampels from float to double precision
 * @param voiceBuffer : buffer with voice signal (of size samplesPerBlock)
 * @param synthBuffer : buffer with synth signal (of size samplesPerBlock)
 */
void MyBuffer::fillInputBuffers(const AudioBuffer<float> &voiceBuffer, const AudioBuffer<float> &synthBuffer)
{
    // For voice buffer
    int numSamples = voiceBuffer.getNumSamples();
    for(int channel = 0; channel < numChannelsInVoice; channel++) {
        auto *voiceReadPtr = voiceBuffer.getReadPointer(channel);
        auto *inVoiceWrtPtr = mInputVoice.getWritePointer(channel, 0);

        #pragma clang loop vectorize(enable)
        for (int i = 0; i < numSamples; i++)
            inVoiceWrtPtr[(inCounter + i) % inSize] = voiceReadPtr[i];

    }

    // For synth buffer
    for(int channel = 0; channel < numChannelsInSynth; channel++) {
        auto* synthReadPtr = synthBuffer.getReadPointer(channel);
        auto* inSynthWrtPtr = mInputSynth.getWritePointer(channel, 0);

        if (synthReadPtr!=NULL) {
            #pragma clang loop vectorize(enable)
            for (int i = 0; i < numSamples; i++)
                inSynthWrtPtr[(inCounter + i) % inSize] = synthReadPtr[i];
        }
        else {
            #pragma clang loop vectorize(enable)
            for (int i = 0; i < numSamples; i++)
                inSynthWrtPtr[(inCounter + i) % inSize] = 0.0;
        }

    }
}


/**
 * Fill the float buffer with content mOutputBuffer. Then updates the counters
 * @param buffer : buffer to copy data to
 * @param nOutputChannels : number of output channels of buffer
 */
void MyBuffer::fillOutputBuffer(AudioBuffer<float> &buffer, int nOutputChannels)
{
    buffer.clear();

    for(int channel = 0; channel < nOutputChannels; channel++) {
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


/**
 * Returns voice sample (double) for given channel and idx
 * @param channel
 * @param idx
 * @return requested samples (double)
 */
double MyBuffer::getVoiceSample(int channel, int idx) const
{
    /*
    if (idx < -samplesToKeep || idx >= (samplesPerBlock + latency)) {
        std::cerr << "idx out of range in getVoiceSample()." << std::endl;
        std::cerr << "idx: " << idx << std::endl;
        assert(false);
    }
     */

    return mInputVoice.getSample(channel, (currCounter + idx + inSize)%inSize);
}


/**
 * Returns synth sample (double) for given channel and index
 * @param channel
 * @param idx
 * @return requested samples (double)
 */
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


/**
 * add a value to the output buffer at given channel and index
 * @param channel
 * @param idx
 * @param value
 */
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


/**
 * returns output sample for a given channel and index
 * @param channel
 * @param idx
 * @return
 */
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


/**
 * Clears numSamples of the given channel of output buffer starting at outCounter
 * @param channel
 * @param numSamples
 */
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

/**
 * Return RMS value of a range of the voice buffer
 * the range start at startSample and is numSamples long
 * @param startSample
 * @param numSamples
 * @return
 */
double MyBuffer::getRMSLevelVoice(int startSample, int numSamples)
{
    int startIdx = (currCounter + startSample + inSize)%inSize;
    double RMS;
    if (startIdx + numSamples <= inSize)
    {
        RMS = mInputVoice.getRMSLevel(0, startIdx, numSamples);
    } else
    {
        RMS = sqrt(pow(mInputVoice.getRMSLevel(0, startIdx, inSize-startIdx), 2) +
                pow(mInputVoice.getRMSLevel(0, 0,
                        numSamples - (inSize - startIdx)), 2));
    }

    return RMS;
}

/**
 * Get RMS level of full voice buffer
 * @return RMS level (double)
 */
double MyBuffer::getRMSLevelVoiceFull()
{
    return mInputVoice.getRMSLevel(0, 0, inSize);
}

/**
 * Return RMS value of a range of the synthesiser buffer
 * the range start at startSample and is numSamples long
 * @param startSample
 * @param numSamples
 * @return
 */
double MyBuffer::getRMSLevelSynth(int startSample, int numSamples)
{
    int startIdx = (currCounter + startSample + inSize)%inSize;
    double RMS_0;
    double RMS_1;
    if (startIdx + numSamples <= inSize)
    {
        RMS_0 = mInputSynth.getRMSLevel(0, startIdx, numSamples);
        RMS_1 = mInputSynth.getRMSLevel(1, startIdx, numSamples);

    } else
    {
        RMS_0 = sqrt(pow(mInputSynth.getRMSLevel(0, startIdx, inSize-startIdx), 2) +
                   pow(mInputSynth.getRMSLevel(0, 0,
                                               numSamples - (inSize - startIdx)), 2));

        RMS_1 = sqrt(pow(mInputSynth.getRMSLevel(1, startIdx, inSize-startIdx), 2) +
                     pow(mInputSynth.getRMSLevel(1, 0,
                                                 numSamples - (inSize - startIdx)), 2));
    }

    // Return mean of RMS value of both channels
    return (RMS_0+RMS_1)/2.0;
}

/**
 * Get RMS level of full synth buffer
 * @return RMS level (double)
 */
double MyBuffer::getRMSLevelSynthFull()
{
    return mInputSynth.getRMSLevel(0, 0, inSize);
}


/**
 * Add dry voice signal with given gain (not in db) to output buffer
 * @param gain: value larger than 0
 */
void MyBuffer::addDryVoice(double gain)
{
    if (currCounter + samplesPerBlock <= inSize)
    {
        if (outCounter + samplesPerBlock <= outSize) 
        {
            for (int channel = 0; channel < numChannelsOut; channel++)
                mOutput.addFrom(channel, outCounter, mInputVoice, 0, currCounter, samplesPerBlock, gain);
        }
        else
        {
            int n = outSize - outCounter;
            for (int channel = 0; channel < numChannelsOut; channel++) 
            {
                mOutput.addFrom(channel, outCounter, mInputVoice, 0, currCounter, n, gain);
                mOutput.addFrom(channel, 0, mInputVoice, 0, currCounter + n,
                        samplesPerBlock - n, gain);
            }

        }
    } 
    else {
        if (outCounter + samplesPerBlock <= outSize) {
            for (int channel = 0; channel < numChannelsOut; channel++) {
                int n = inSize - currCounter;
                mOutput.addFrom(channel, outCounter, mInputVoice, 0, currCounter, n, gain);
                mOutput.addFrom(channel, outCounter + n, mInputVoice, 0, 0,
                                samplesPerBlock - n, gain);
            }
        } else {
            int nOut = outSize - outCounter;
            int nIn =  inSize - currCounter;

            if (nIn == nOut) {
                for (int channel = 0; channel < numChannelsOut; channel++) {
                    mOutput.addFrom(channel, outCounter, mInputVoice, 0, currCounter, nIn, gain);
                    mOutput.addFrom(channel, 0, mInputVoice, 0, 0,
                            samplesPerBlock - nIn, gain);
                }
            }

            if (nIn < nOut) {
                for (int channel = 0; channel < numChannelsOut; channel++) {
                    mOutput.addFrom(channel, outCounter, mInputVoice, 0, currCounter, nIn, gain);
                    mOutput.addFrom(channel, outCounter + nIn, mInputVoice, 0, 
                            0, nOut - nIn, gain);
                    mOutput.addFrom(channel, 0, mInputVoice, 0, nOut - nIn,
                                    samplesPerBlock - nOut, gain);

                }
            }

            if (nIn > nOut) {
                for (int channel = 0; channel < numChannelsOut; channel++) {
                    mOutput.addFrom(channel, outCounter, mInputVoice, 0, currCounter, nOut, gain);
                    mOutput.addFrom(channel, 0, mInputVoice, 0, currCounter + nOut,
                                    nIn - nOut, gain);
                    mOutput.addFrom(channel, nIn - nOut, mInputVoice, 0, 0,
                                    samplesPerBlock - nIn, gain);
                }
            }
        }
    } 

}


/**
 * Add synth signal with given gain (not in db) to output buffer
 * @param gain: value larger than 0
 */
void MyBuffer::addSynth(double gain)
{
    if (currCounter + samplesPerBlock <= inSize)
    {
        if (outCounter + samplesPerBlock <= outSize)
        {
            for (int channel = 0; channel < numChannelsOut; channel++)
                mOutput.addFrom(channel, outCounter, mInputSynth, channel, currCounter, samplesPerBlock, gain);
        }
        else
        {
            int n = outSize - outCounter;
            for (int channel = 0; channel < numChannelsOut; channel++)
            {
                mOutput.addFrom(channel, outCounter, mInputSynth, channel, currCounter, n, gain);
                mOutput.addFrom(channel, 0, mInputSynth, channel, currCounter + n,
                        samplesPerBlock - n, gain);
            }

        }
    }
    else {
        if (outCounter + samplesPerBlock <= outSize) {
            for (int channel = 0; channel < numChannelsOut; channel++)
            {
                int n = inSize - currCounter;
                mOutput.addFrom(channel, outCounter, mInputSynth, channel, currCounter, n, gain);
                mOutput.addFrom(channel, outCounter + n, mInputSynth, channel, 0,
                                samplesPerBlock - n, gain);
            }
        } else {
            int nOut = outSize - outCounter;
            int nIn =  inSize - currCounter;

            if (nIn == nOut) {
                for (int channel = 0; channel < numChannelsOut; channel++)
                {
                    mOutput.addFrom(channel, outCounter, mInputSynth, channel, currCounter, nIn, gain);
                    mOutput.addFrom(channel, 0, mInputSynth, channel, 0,
                            samplesPerBlock - nIn, gain);
                }
            }

            if (nIn < nOut) {
                for (int channel = 0; channel < numChannelsOut; channel++)
                {
                    mOutput.addFrom(channel, outCounter, mInputSynth, channel, currCounter, nIn, gain);
                    mOutput.addFrom(channel, outCounter + nIn, mInputSynth, channel,
                                    0, nOut - nIn, gain);
                    mOutput.addFrom(channel, 0, mInputSynth, channel, nOut - nIn,
                                    samplesPerBlock - nOut, gain);

                }
            }

            if (nIn > nOut) {
                for (int channel = 0; channel < numChannelsOut; channel++)
                {
                    mOutput.addFrom(channel, outCounter, mInputSynth, channel, currCounter, nOut, gain);
                    mOutput.addFrom(channel, 0, mInputSynth, channel, currCounter + nOut,
                                    nIn - nOut, gain);
                    mOutput.addFrom(channel, nIn - nOut, mInputSynth, channel, 0,
                                    samplesPerBlock - nIn, gain);
                }
            }
        }
    }

}
