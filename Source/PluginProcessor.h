/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include "../JuceLibraryCode/JuceHeader.h"
#include "MyBuffer.h"
#include "VocoderProcess.h"
#include "PitchProcess.h"
#include <string>
#include <math.h>


//==============================================================================
/**
*/
class VocoderAudioProcessor  : public AudioProcessor
{
public:
    //==============================================================================
    VocoderAudioProcessor();
    ~VocoderAudioProcessor();

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

    void processBlock (AudioBuffer<float>&, MidiBuffer&) override;

    //==============================================================================
    AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const String getProgramName (int index) override;
    void changeProgramName (int index, const String& newName) override;

    //==============================================================================
    void getStateInformation (MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

    AudioProcessorValueTreeState::ParameterLayout createParameterLayout();
    AudioProcessorValueTreeState treeState;

private:
    // To handle the variable buffer size and latency
    MyBuffer myBuffer;

    // Class with process function
    VocoderProcess vocoderProcess;
    PitchProcess pitchProcess;

    // GUI
    foleys::MagicProcessorState magicState { *this, treeState };

    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (VocoderAudioProcessor)
};
