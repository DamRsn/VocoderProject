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
                       .withInput  ("Input",  AudioChannelSet::mono(), true)

#endif
                       .withOutput ("Output", AudioChannelSet::stereo(), true)
                       .withInput  ("Sidechain",  AudioChannelSet::stereo())
#endif
                       ),
treeState(*this, nullptr, "PARAMETERS", createParameterLayout())
#endif
{
    vocoderProcess.setAudioProcPtr(this);
    pitchProcess.setAudioProcPtr(this);
    //orderMaxVoice = 100;
    //orderMaxSynth = 30;

}

AudioProcessorValueTreeState::ParameterLayout VocoderAudioProcessor::createParameterLayout()
{
    std::vector <std::unique_ptr<RangedAudioParameter>> params;

    auto gainPitchParam = std::make_unique<AudioParameterFloat>("gainPitch", "GainPitch", -60.0f, 6.0f, 0.0f);
    params.push_back(std::move(gainPitchParam));

    auto gainVoiceParam = std::make_unique<AudioParameterFloat>("gainVoice", "GainVoice", -60.0f, 6.0f, 0.0f);
    params.push_back(std::move(gainVoiceParam));

    auto gainSynthParam = std::make_unique<AudioParameterFloat>("gainSynth", "GainSynth", -60.0f, 6.0f, 0.0f);
    params.push_back(std::move(gainSynthParam));

    auto gainVocParam = std::make_unique<AudioParameterFloat>("gainVoc", "GainVoc", -60.0f, 6.0f, 0.0f);
    params.push_back(std::move(gainVocParam));

    auto lpcVoiceParam = std::make_unique<AudioParameterInt>("lpcVoice", "LpcVoice", 2, 100, 30);
    params.push_back(std::move(lpcVoiceParam));

    auto lpcSynthParam = std::make_unique<AudioParameterInt>("lpcSynth", "LpcSynth", 2, 30, 5);
    params.push_back(std::move(lpcSynthParam));

    return { params.begin(), params.end() };
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
    /*
     gainSynth = -60.0;
     gainVocoder = -60.0;
     gainVoice = -60.0;
     orderVoice = 20;
     orderSynth = 5;
     orderMaxVoice = 100;
     orderMaxSynth = 20;
    */

    auto totalNumInputChannels = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();

    if (totalNumInputChannels != 3 && totalNumOutputChannels != 2)
    {
        std::cerr << "Wrong number of input or output channel. " << std::endl;
        assert(false);
    }

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

    /*
    orderVoice = 15;
    orderSynth = 4;
     */

    pitchProcess.prepare(sampleRate, 100, 800, frameLenPitch, hopPitch, 10, 40,
                         0, samplesPerBlock, Notes::Chrom);
    vocoderProcess.prepare(wlenVoc, hopVoc, 0, 0, window_str);

    int latency = std::max(pitchProcess.getLatency(samplesPerBlock), vocoderProcess.getLatency(samplesPerBlock));
    int samplesToKeep = frameLenPitch;

    // TODO: this was here just for test
    //latency = vocoderProcess.getLatency(samplesPerBlock);
    //samplesToKeep = 0;


    myBuffer.prepare(samplesPerBlock, samplesToKeep, latency, sampleRate, 1, 2,
            totalNumOutputChannels);

    pitchProcess.prepare2(myBuffer);

    setLatencySamples(latency);
}

void VocoderAudioProcessor::releaseResources() {}

#ifndef JucePlugin_PreferredChannelConfigurations
bool VocoderAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
    if (layouts.getMainInputChannelSet() == AudioChannelSet::mono() &&
    layouts.getMainOutputChannelSet() == AudioChannelSet::stereo() &&
    !layouts.getChannelSet(true, 1).isDisabled())
    {
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
