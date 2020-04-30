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
    gainSlider.setTextBoxStyle(Slider::TextBoxBelow, true, 100, 25);
    gainSlider.setRange(0.0f, 5.0f, 0.01f);
    gainSlider.setValue(1.0);
    gainSlider.addListener(this);
    addAndMakeVisible(gainSlider);

    velocitySlider.setSliderStyle(Slider::SliderStyle::LinearVertical);
    velocitySlider.setTextBoxStyle(Slider::TextBoxBelow, true, 100, 25);
    velocitySlider.setRange(0.0f, 1.0f, 0.01f);
    velocitySlider.setValue(0.2);
    velocitySlider.addListener(this);
    addAndMakeVisible(velocitySlider);

    setSize (400, 300);

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
    gainSlider.setBoundsRelative(0.25f, 0.0f, 1.0f, 1.0f);
    velocitySlider.setBoundsRelative(-0.25f, 0.0f,1.0f, 1.0f);

}


void VocoderAudioProcessorEditor::sliderValueChanged(Slider *slider)
{
    if (slider == &gainSlider)
    {
        processor.gain = gainSlider.getValue();
    }

    if (slider == &velocitySlider)
    {
        processor.velocity = velocitySlider.getValue();
    }
}