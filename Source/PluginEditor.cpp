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
    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.

    gainSlider.setSliderStyle(Slider::SliderStyle::LinearVertical);
    gainSlider.setTextBoxStyle(Slider::TextBoxBelow, true, 50, 25);
    gainSlider.setRange(-60.0f, 60.0f, 0.1f);
    gainSlider.setValue(0);
    gainSlider.addListener(this);
    addAndMakeVisible(gainSlider);
    gainLabel.setText("G bIIR", dontSendNotification);
    gainLabel.attachToComponent(&gainSlider, true);
    addAndMakeVisible(gainLabel);

    /*
    velocitySlider.setSliderStyle(Slider::SliderStyle::LinearVertical);
    velocitySlider.setTextBoxStyle(Slider::TextBoxBelow, true, 50, 25);
    velocitySlider.setRange(0.0f, 1.0f, 0.01f);
    velocitySlider.setValue(0.0);
    velocitySlider.addListener(this);
    addAndMakeVisible(velocitySlider);
    velocityLabel.setText("Vel NW", dontSendNotification);
    velocityLabel.attachToComponent(&velocitySlider, true);
    addAndMakeVisible(velocityLabel);
    */
    gainVoiceSlider.setSliderStyle(Slider::SliderStyle::LinearVertical);
    gainVoiceSlider.setTextBoxStyle(Slider::TextBoxBelow, true, 50, 25);
    gainVoiceSlider.setRange(-60.0f, 60.0f, 0.1f);
    gainVoiceSlider.setValue(-60);
    gainVoiceSlider.addListener(this);
    addAndMakeVisible(gainVoiceSlider);
    gainVoiceLabel.setText("G Voice", dontSendNotification);
    gainVoiceLabel.attachToComponent(&gainVoiceSlider, true);
    addAndMakeVisible(gainVoiceLabel);
    
    gainSynthSlider.setSliderStyle(Slider::SliderStyle::LinearVertical);
    gainSynthSlider.setTextBoxStyle(Slider::TextBoxBelow, true, 50, 25);
    gainSynthSlider.setRange(-60.0f, 60.0f, 0.1f);
    gainSynthSlider.setValue(-60);
    gainSynthSlider.addListener(this);
    addAndMakeVisible(gainSynthSlider);
    gainSynthLabel.setText("G Synth", dontSendNotification);
    gainSynthLabel.attachToComponent(&gainSynthSlider, true);
    addAndMakeVisible(gainSynthLabel);
    
    gainVocoderSlider.setSliderStyle(Slider::SliderStyle::LinearVertical);
    gainVocoderSlider.setTextBoxStyle(Slider::TextBoxBelow, true, 50, 25);
    gainVocoderSlider.setRange(-60.0f, 60.0f, 0.1f);
    gainVocoderSlider.setValue(0.0);
    gainVocoderSlider.addListener(this);
    addAndMakeVisible(gainVocoderSlider);
    gainVocoderLabel.setText("G Voc", dontSendNotification);
    gainVocoderLabel.attachToComponent(&gainVocoderSlider, true);
    addAndMakeVisible(gainVocoderLabel);


    LPCVoiceSlider.setSliderStyle(Slider::SliderStyle::LinearVertical);
    LPCVoiceSlider.setTextBoxStyle(Slider::TextBoxBelow, true, 50, 25);
    LPCVoiceSlider.setRange(2, processor.orderMaxVoice, 1);
    LPCVoiceSlider.setValue(15);
    LPCVoiceSlider.addListener(this);
    addAndMakeVisible(LPCVoiceSlider);
    LPCVoiceLabel.setText("LPC Vx", dontSendNotification);
    LPCVoiceLabel.attachToComponent(&LPCVoiceSlider, true);
    addAndMakeVisible(LPCVoiceLabel);

    LPCSynthSlider.setSliderStyle(Slider::SliderStyle::LinearVertical);
    LPCSynthSlider.setTextBoxStyle(Slider::TextBoxBelow, true, 50, 25);
    LPCSynthSlider.setRange(2, processor.orderMaxSynth, 1);
    LPCSynthSlider.setValue(4);
    LPCSynthSlider.addListener(this);
    addAndMakeVisible(LPCSynthSlider);
    LPCSynthLabel.setText("LPC Sth", dontSendNotification);
    LPCSynthLabel.attachToComponent(&LPCSynthSlider, true);
    addAndMakeVisible(LPCSynthLabel);
    
    
    

    setSize (1000, 300);

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
    gainSlider.setBoundsRelative(0.1f, 0.0f, 0.1f, 0.8f);
    //velocitySlider.setBoundsRelative(0.23f, 0.0f,0.1f, 0.8f);
    gainVoiceSlider.setBoundsRelative(0.26f, 0.0f, 0.1f, 0.8f);
    gainSynthSlider.setBoundsRelative(0.42f, 0.0f, 0.1f, 0.8f);
    gainVocoderSlider.setBoundsRelative(0.58f, 0.0f, 0.1f, 0.8f);
    LPCSynthSlider.setBoundsRelative(0.74f, 0.0f, 0.1f, 0.8f);
    LPCVoiceSlider.setBoundsRelative(0.9f, 0.0f, 0.1f, 0.8f);

}


void VocoderAudioProcessorEditor::sliderValueChanged(Slider *slider)
{
    if (slider == &gainSlider)
        processor.gain = gainSlider.getValue();

    //if (slider == &velocitySlider)
    //    processor.velocity = velocitySlider.getValue();

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