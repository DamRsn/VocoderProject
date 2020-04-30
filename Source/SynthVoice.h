/*
  ==============================================================================

    SynthVoice.h
    Created: 15 Apr 2020 11:52:43am
    Author:  Damien Ronssin

  ==============================================================================
*/

#pragma once
#include "../JuceLibraryCode/JuceHeader.h"
#include "SynthSound.h"
#include "maximilian.h"
#include <iostream>

class VocoderAudioProcessor;


class SynthVoice : public SynthesiserVoice
{
public:

    //SynthVoice(VocoderAudioProcessor* audioProcPtr);
    bool canPlaySound (SynthesiserSound* sound) override;

    void startNote (int midiNoteNumber, float velocity, SynthesiserSound *sound, int currentPitchWheelPosition)
    override;
    
    void stopNote (float velocity, bool allowTailOff) override;
    
    void pitchWheelMoved (int newPitchWheelValue) override;
    
    void controllerMoved (int controllerNumber, int newControllerValue) override;
    
    void renderNextBlock (AudioBuffer< float > &outputBuffer, int startSample, int numSamples) override;

private:
    double level;
    double frequency;
    
    maxiOsc  osc1;
    maxiEnv env1;
    //VocoderAudioProcessor* audioProcPtr;
};

