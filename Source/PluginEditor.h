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
class VocoderAudioProcessorEditor  : public AudioProcessorEditor
{
public:
    VocoderAudioProcessorEditor (VocoderAudioProcessor&);
    ~VocoderAudioProcessorEditor();

    //==============================================================================
    void paint (Graphics&) override;
    void resized() override;

private:
    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    VocoderAudioProcessor& processor;

    Slider gainPitchSlider;
    Label gainPitchLabel;

    Slider gainVoiceSlider;
    Label gainVoiceLabel;

    Slider gainSynthSlider;
    Label gainSynthLabel;

    Slider gainVocoderSlider;
    Label gainVocoderLabel;

    Slider LPCVoiceSlider;
    Label LPCVoiceLabel;

    Slider LPCSynthSlider;
    Label LPCSynthLabel;

    ComboBox keyBox;
    Label keyBoxLabel;

    TextButton pitchButton;

    TextButton vocButton;


// another public for destruction order
public:

    std::unique_ptr<AudioProcessorValueTreeState::SliderAttachment> pitchSliderValue;
    std::unique_ptr<AudioProcessorValueTreeState::SliderAttachment> voiceSliderValue;
    std::unique_ptr<AudioProcessorValueTreeState::SliderAttachment> synthSliderValue;
    std::unique_ptr<AudioProcessorValueTreeState::SliderAttachment> vocoderSliderValue;
    std::unique_ptr<AudioProcessorValueTreeState::SliderAttachment> LPCVoiceSliderValue;
    std::unique_ptr<AudioProcessorValueTreeState::SliderAttachment> LPCSynthSliderValue;

    std::unique_ptr<AudioProcessorValueTreeState::ComboBoxAttachment> keyBoxValue;
    std::unique_ptr<AudioProcessorValueTreeState::ButtonAttachment> pitchButtonValue;
    std::unique_ptr<AudioProcessorValueTreeState::ButtonAttachment> vocButtonValue;




    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (VocoderAudioProcessorEditor)
};
