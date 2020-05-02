/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include "../JuceLibraryCode/JuceHeader.h"
#include "PluginProcessor.h"

//==============================================================================
/**
*/
class VocoderAudioProcessorEditor  : public AudioProcessorEditor, public Slider::Listener
{
public:
    VocoderAudioProcessorEditor (VocoderAudioProcessor&);
    ~VocoderAudioProcessorEditor();

    //==============================================================================
    void paint (Graphics&) override;
    void resized() override;

    void sliderValueChanged(Slider* slider) override;

private:
    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    VocoderAudioProcessor& processor;

    Slider gainSlider;
    Slider velocitySlider;
    Slider gainVoiceSlider;
    Slider gainSynthSlider;
    Slider gainVocoderSlider;
    Slider LPCVoiceSlider;
    Slider LPCSynthSlider;

    Label gainLabel;
    Label velocityLabel;
    Label gainVoiceLabel;
    Label gainSynthLabel;
    Label gainVocoderLabel;
    Label LPCVoiceLabel;
    Label LPCSynthLabel;



    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (VocoderAudioProcessorEditor)
};
