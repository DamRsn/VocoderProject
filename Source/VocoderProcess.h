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


class VocoderProcess
{
public:
    VocoderProcess();
    void prepare(int wlen, int hop, std::string windowType);
    ~VocoderProcess();

    void process(MyBuffer& myBuffer);
    std::vector<float> window;


private:
    int wlen;
    int hop;
    int startSample;

    int k;

    void processWindow(MyBuffer& myBuffer);




};