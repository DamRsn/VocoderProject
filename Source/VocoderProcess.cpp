/*
  ==============================================================================

    VocoderProcess.cpp
    Created: 20 Apr 2020 10:17:42am
    Author:  Damien Ronssin

  ==============================================================================
*/

#include "VocoderProcess.h"

VocoderProcess::VocoderProcess() {    
}

void VocoderProcess::prepare(int wlen_, int hop_, std::string windowType) {

    wlen = wlen_;
    hop = hop_;
    startSample = 0;
    window.resize(wlen, 0.0);
    k = 0;
    
    if (windowType=="hann")
        dsp::WindowingFunction<float>::fillWindowingTables(&window[0], wlen, dsp::WindowingFunction<float>::hann, false);
    else
        std::cerr << "Unknown window type" << std::endl;

}

VocoderProcess::~VocoderProcess() {}


void VocoderProcess::process(MyBuffer &myBuffer)
{
    /*
    if (k < 5)
        std::cout<< "SS before while: " << startSample << std::endl;
     */
    while (startSample < myBuffer.getSamplesPerBlock())
    {
        /*
        if (k < 5)
            std::cout<< "SS: " << startSample << std::endl;
         */
        processWindow(myBuffer);
        startSample += hop;
    }

    startSample -= myBuffer.getSamplesPerBlock();
    k +=1;
    k %= 5000;
}

void VocoderProcess::processWindow(MyBuffer &myBuffer)
{
    int numChannels = myBuffer.getNumChannels();
    float value;

    for (int channel = 0; channel < numChannels; channel++)
    {
        for (int i = startSample; i < (startSample + wlen); i++)
        {
            //value = (myBuffer.getSynthSample(channel, i) + myBuffer.getVoiceSample(channel, i)) * window[i-startSample];
            value = myBuffer.getVoiceSample(channel, i) * window[i-startSample];
            myBuffer.addOutSample(channel, i, value);
        }
    }
}

