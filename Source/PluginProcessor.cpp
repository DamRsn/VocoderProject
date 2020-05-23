/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
VocoderAudioProcessor::VocoderAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
     : AudioProcessor (BusesProperties()
                     #if ! JucePlugin_IsMidiEffect
                      #if ! JucePlugin_IsSynth
                       .withInput  ("Input",  AudioChannelSet::stereo(), true)

#endif
                       .withOutput ("Output", AudioChannelSet::stereo(), true)
                       .withInput  ("Sidechain",  AudioChannelSet::stereo())


#endif
                       )
#endif
{

    vocoderProcess.setAudioProcPtr(this);
    pitchProcess.setAudioProcPtr(this);
    mySynth.clearVoices();
    for (int i = 0; i < 10; i++)
    {
        mySynth.addVoice(new SynthVoice());
    }
    mySynth.clearSounds();
    
    mySynth.addSound(new SynthSound());
    orderMaxVoice = 100;
    orderMaxSynth = 30;

}

VocoderAudioProcessor::~VocoderAudioProcessor()
{
}

//==============================================================================
const String VocoderAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool VocoderAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool VocoderAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool VocoderAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double VocoderAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int VocoderAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int VocoderAudioProcessor::getCurrentProgram()
{
    return 0;
}

void VocoderAudioProcessor::setCurrentProgram (int index)
{
}

const String VocoderAudioProcessor::getProgramName (int index)
{
    return {};
}

void VocoderAudioProcessor::changeProgramName (int index, const String& newName)
{
}

//==============================================================================
void VocoderAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock) {
    ignoreUnused(samplesPerBlock);
    lastSampleRate = sampleRate;
    mySynth.setCurrentPlaybackSampleRate(lastSampleRate);

    synthBuffer.setSize(2, samplesPerBlock);

    auto totalNumInputChannels = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();

    // Ratio of sampling rate to get correct size for different arrays ...
    double ratioSR = sampleRate/44100.0;

    // Vocoder
    int hopVoc = floor(128.0 * ratioSR);
    int wlenVoc = 4*hopVoc;
    std::string window_str = "sine";
    int orderMaxVoc = 100;


    // Pitch Corrector
    int corres_256 = floor(256.0 * ratioSR);
    int hopPitch = 3 * corres_256;
    int frameLenPitch = 4 * corres_256;

    orderVoice = 15;
    orderSynth = 4;

    pitchProcess.prepare(sampleRate, 100, 800, frameLenPitch, hopPitch, 10, 15,
                         0, samplesPerBlock, Notes::Chrom);
    vocoderProcess.prepare(wlenVoc, hopVoc, orderVoice, orderMaxVoice, window_str);

    int latency = std::max(pitchProcess.getLatency(samplesPerBlock), vocoderProcess.getLatency(samplesPerBlock));
    int samplesToKeep = floor(1024.0 * ratioSR);

    myBuffer.prepare(samplesPerBlock, samplesToKeep, latency, sampleRate, totalNumOutputChannels);

    setLatencySamples(latency);

}

void VocoderAudioProcessor::releaseResources() {}

#ifndef JucePlugin_PreferredChannelConfigurations
bool VocoderAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
    if (layouts.getMainInputChannelSet() == layouts.getMainOutputChannelSet() && ! layouts.getChannelSet(true, 1).isDisabled()) {
        return true;
    }
    else {
        return false;
    }
}
#endif

void VocoderAudioProcessor::processBlock (AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
{
    ScopedNoDenormals noDenormals;

    auto samplesPerBlock = buffer.getNumSamples();
    auto nOutputChannels = getTotalNumOutputChannels();
    auto voiceBuffer = getBusBuffer(buffer, true, 0);
    auto synthBuffer = getBusBuffer(buffer, true, 1);

    /*
    for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
        buffer.clear (i, 0, buffer.getNumSamples());

    synthBuffer.clear();
    mySynth.renderNextBlock(synthBuffer, midiMessages, 0, synthBuffer.getNumSamples());
    */

    myBuffer.fillInputBuffers(voiceBuffer, synthBuffer);

    vocoderProcess.process(myBuffer);
    pitchProcess.process(myBuffer);

    myBuffer.fillOutputBuffer(buffer, nOutputChannels);
}

//==============================================================================
bool VocoderAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

AudioProcessorEditor* VocoderAudioProcessor::createEditor()
{
    return new VocoderAudioProcessorEditor (*this);
}

//==============================================================================
void VocoderAudioProcessor::getStateInformation (MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.
}

void VocoderAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.
}

//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new VocoderAudioProcessor();
}
