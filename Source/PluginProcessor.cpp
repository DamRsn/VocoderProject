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
                     #endif
                       )
#endif
{
    mySynth.clearVoices();
    for (int i = 0; i < 10; i++)
    {
        mySynth.addVoice(new SynthVoice());
    }
    mySynth.clearSounds();
    
    mySynth.addSound(new SynthSound());
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
void VocoderAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    ignoreUnused(samplesPerBlock);
    lastSampleRate = sampleRate;
    mySynth.setCurrentPlaybackSampleRate(lastSampleRate);

    synthBuffer.setSize(2, samplesPerBlock);

    auto totalNumInputChannels  = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();

    int wlen = 256;
    int hop = 128;
    std::string window_str = "hann";
    
    int samplesToKeep = 0;
    int latency;
    
    if(samplesPerBlock<= hop)
        latency = wlen - samplesPerBlock;
    else
        latency = wlen;
    
    // latency = int(ceil(128.0/double(samplesPerBlock))) * samplesPerBlock;

    std::cout << "\nbuffersize:  " << samplesPerBlock << std::endl;
    std::cout << "inSize:  " << samplesPerBlock + samplesToKeep + latency << std::endl;
    std::cout << "samplesToKeep:  " << samplesToKeep << std::endl;
    std::cout << "latency:  " << latency << std::endl;
    std::cout << "wlen:  " << wlen << std::endl;
    std::cout << "hop:  " << hop << std::endl;

    

    vocoderProcess.prepare(wlen, hop, window_str);

    myBuffer.prepare(samplesPerBlock, samplesToKeep, latency, sampleRate, getTotalNumInputChannels());

    setLatencySamples(latency);

    k = 0;

    /*
    for (int i=0; i < wlen; i++)
    {
        std::cout<< vocoderProcess.window[i] << ", ";
    }
    */
}

void VocoderAudioProcessor::releaseResources()
{
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool VocoderAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
  #if JucePlugin_IsMidiEffect
    ignoreUnused (layouts);
    return true;
  #else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    if (layouts.getMainOutputChannelSet() != AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
   #if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
   #endif

    return true;
  #endif
}
#endif

void VocoderAudioProcessor::processBlock (AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
{
    ScopedNoDenormals noDenormals;
    
    auto totalNumInputChannels  = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();
    auto samplesPerBlock = buffer.getNumSamples();
    
    for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
        buffer.clear (i, 0, buffer.getNumSamples());

    synthBuffer.clear();
    mySynth.renderNextBlock(synthBuffer, midiMessages, 0, synthBuffer.getNumSamples());

    myBuffer.fillInputBuffers(buffer, synthBuffer);

    vocoderProcess.process(myBuffer);

    myBuffer.fillOutputBuffer(buffer);
    
    
    
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
