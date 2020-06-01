/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
VocoderAudioProcessorEditor::VocoderAudioProcessorEditor (VocoderAudioProcessor& p)
    : AudioProcessorEditor (&p), processor (p)
{
    pitchSliderValue = std::make_unique<AudioProcessorValueTreeState::SliderAttachment>(processor.treeState,
            "gainPitch", gainPitchSlider);

    voiceSliderValue = std::make_unique<AudioProcessorValueTreeState::SliderAttachment>(processor.treeState,
                                                                                        "gainVoice", gainVoiceSlider);

    vocoderSliderValue = std::make_unique<AudioProcessorValueTreeState::SliderAttachment>(processor.treeState,
                                                                                        "gainVoc", gainVocoderSlider);

    synthSliderValue = std::make_unique<AudioProcessorValueTreeState::SliderAttachment>(processor.treeState,
                                                                                        "gainSynth", gainSynthSlider);

    LPCVoiceSliderValue = std::make_unique<AudioProcessorValueTreeState::SliderAttachment>(processor.treeState,
                                                                                        "lpcVoice", LPCVoiceSlider);

    LPCSynthSliderValue = std::make_unique<AudioProcessorValueTreeState::SliderAttachment>(processor.treeState,
                                                                                        "lpcSynth", LPCSynthSlider);

    keyBoxValue = std::make_unique<AudioProcessorValueTreeState::ComboBoxAttachment>(processor.treeState,
                                                                                           "keyPitch", keyBox);

    pitchButtonValue = std::make_unique<AudioProcessorValueTreeState::ButtonAttachment>(processor.treeState,
                                                                                     "pitchBool", pitchButton);

    vocButtonValue = std::make_unique<AudioProcessorValueTreeState::ButtonAttachment>(processor.treeState,
                                                                                        "vocBool", vocButton);



    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.

    gainPitchSlider.setSliderStyle(Slider::SliderStyle::LinearVertical);
    gainPitchSlider.setTextBoxStyle(Slider::TextBoxBelow, true, 50, 25);
    //gainPitchSlider.setRange(-60.0f, 6.0f, 0.1f);
    //gainPitchSlider.setValue(0);
    //gainPitchSlider.addListener(this);
    addAndMakeVisible(gainPitchSlider);
    gainPitchLabel.setText("G Pitch", dontSendNotification);
    gainPitchLabel.attachToComponent(&gainPitchSlider, true);
    addAndMakeVisible(gainPitchLabel);


    gainVoiceSlider.setSliderStyle(Slider::SliderStyle::LinearVertical);
    gainVoiceSlider.setTextBoxStyle(Slider::TextBoxBelow, true, 50, 25);
    //gainVoiceSlider.setRange(-60.0f, 6.0f, 0.1f);
    //gainVoiceSlider.setValue(-60);
    //gainVoiceSlider.addListener(this);
    addAndMakeVisible(gainVoiceSlider);
    gainVoiceLabel.setText("G Voice", dontSendNotification);
    gainVoiceLabel.attachToComponent(&gainVoiceSlider, true);
    addAndMakeVisible(gainVoiceLabel);
    
    gainSynthSlider.setSliderStyle(Slider::SliderStyle::LinearVertical);
    gainSynthSlider.setTextBoxStyle(Slider::TextBoxBelow, true, 50, 25);
    //gainSynthSlider.setRange(-60.0f, 6.0f, 0.1f);
    //gainSynthSlider.setValue(-60);
    //gainSynthSlider.addListener(this);
    addAndMakeVisible(gainSynthSlider);
    gainSynthLabel.setText("G Synth", dontSendNotification);
    gainSynthLabel.attachToComponent(&gainSynthSlider, true);
    addAndMakeVisible(gainSynthLabel);
    
    gainVocoderSlider.setSliderStyle(Slider::SliderStyle::LinearVertical);
    gainVocoderSlider.setTextBoxStyle(Slider::TextBoxBelow, true, 50, 25);
    //gainVocoderSlider.setRange(-60.0f, 6.0f, 0.1f);
    //gainVocoderSlider.setValue(0.0);
    //gainVocoderSlider.addListener(this);
    addAndMakeVisible(gainVocoderSlider);
    gainVocoderLabel.setText("G Voc", dontSendNotification);
    gainVocoderLabel.attachToComponent(&gainVocoderSlider, true);
    addAndMakeVisible(gainVocoderLabel);


    LPCVoiceSlider.setSliderStyle(Slider::SliderStyle::LinearVertical);
    LPCVoiceSlider.setTextBoxStyle(Slider::TextBoxBelow, true, 50, 25);
    //LPCVoiceSlider.setRange(2, processor.orderMaxVoice, 1);
    //LPCVoiceSlider.setValue(15);
    //LPCVoiceSlider.addListener(this);
    addAndMakeVisible(LPCVoiceSlider);
    LPCVoiceLabel.setText("LPC Vx", dontSendNotification);
    LPCVoiceLabel.attachToComponent(&LPCVoiceSlider, true);
    addAndMakeVisible(LPCVoiceLabel);

    LPCSynthSlider.setSliderStyle(Slider::SliderStyle::LinearVertical);
    LPCSynthSlider.setTextBoxStyle(Slider::TextBoxBelow, true, 50, 25);
    //LPCSynthSlider.setRange(2, processor.orderMaxSynth, 1);
    //LPCSynthSlider.setValue(4);
    //LPCSynthSlider.addListener(this);
    addAndMakeVisible(LPCSynthSlider);
    LPCSynthLabel.setText("LPC Sth", dontSendNotification);
    LPCSynthLabel.attachToComponent(&LPCSynthSlider, true);
    addAndMakeVisible(LPCSynthLabel);

    keyBox.addItemList(StringArray("A", "A#", "B", "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "Chrom"),
            0);
    keyBox.setJustificationType (Justification::centred);
    addAndMakeVisible (keyBox);
    keyBoxLabel.setText("Key", dontSendNotification);
    keyBoxLabel.attachToComponent(&keyBox, true);
    addAndMakeVisible(keyBoxLabel);


    pitchButton.setColour (TextButton::buttonColourId, Colours::red);
    pitchButton.setColour (TextButton::buttonOnColourId, Colours::green);
    pitchButton.setClickingTogglesState (true);
    pitchButton.onClick = [this]() { };
    addAndMakeVisible (pitchButton);



    vocButton.setColour (TextButton::buttonColourId, Colours::red);
    vocButton.setColour (TextButton::buttonOnColourId, Colours::green);
    vocButton.setClickingTogglesState (true);

    vocButton.onClick = [this]() { };
    addAndMakeVisible (vocButton);
    

    setSize (900, 300);

}

VocoderAudioProcessorEditor::~VocoderAudioProcessorEditor()
{
}

//==============================================================================
void VocoderAudioProcessorEditor::paint (Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));

    /*
    g.setColour (Colours::white);
    g.setFont (15.0f);
    g.drawFittedText ("Hello World!", getLocalBounds(), Justification::centred, 1);
     */
}

void VocoderAudioProcessorEditor::resized()
{
    gainPitchSlider.setBoundsRelative(0.1f, 0.0f, 0.1f, 0.8f);

    gainVoiceSlider.setBoundsRelative(0.2f, 0.0f, 0.1f, 0.8f);
    gainSynthSlider.setBoundsRelative(0.3f, 0.0f, 0.1f, 0.8f);
    gainVocoderSlider.setBoundsRelative(0.4f, 0.0f, 0.1f, 0.8f);
    LPCSynthSlider.setBoundsRelative(0.5f, 0.0f, 0.1f, 0.8f);
    LPCVoiceSlider.setBoundsRelative(0.6f, 0.0f, 0.1f, 0.8f);

    keyBox.setBoundsRelative(0.7f, 0.0f, 0.1f, 0.8f);

    pitchButton.setBoundsRelative(0.8f, 0.0f, 0.1f, 0.8f);

    vocButton.setBoundsRelative(0.9f, 0.0f, 0.1f, 0.8f);

}

/*
void VocoderAudioProcessorEditor::sliderValueChanged(Slider *slider)
{
    if (slider == &gainPitchSlider)
        processor.gainPitch = gainPitchSlider.getValue();

    if (slider == &gainSynthSlider)
        processor.gainSynth = gainSynthSlider.getValue();

    if (slider == &gainVoiceSlider)
        processor.gainVoice = gainVoiceSlider.getValue();

    if (slider == &gainVocoderSlider)
        processor.gainVocoder = gainVocoderSlider.getValue();

    if (slider == &LPCVoiceSlider)
        processor.orderVoice = LPCVoiceSlider.getValue();

    if (slider == &LPCSynthSlider)
        processor.orderSynth = LPCSynthSlider.getValue();

}
*/
