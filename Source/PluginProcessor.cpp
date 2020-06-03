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
    // Set the pointer to audioProcessor (this) for vocoderProcess and pitchProcess
    vocoderProcess.setAudioProcPtr(this);
    pitchProcess.setAudioProcPtr(this);
}

AudioProcessorValueTreeState::ParameterLayout VocoderAudioProcessor::createParameterLayout()
{
    /*
     * Create all the parameters for the value tree state
     */

    std::vector <std::unique_ptr<RangedAudioParameter>> params;

    auto gainPitchParam = std::make_unique<AudioParameterFloat>("gainPitch", "GainPitch", -60.0f, 6.0f, 0.0f);
    params.push_back(std::move(gainPitchParam));

    auto gainVoiceParam = std::make_unique<AudioParameterFloat>("gainVoice", "GainVoice", -60.0f, 6.0f, -60.0f);
    params.push_back(std::move(gainVoiceParam));

    auto gainSynthParam = std::make_unique<AudioParameterFloat>("gainSynth", "GainSynth", -60.0f, 6.0f, -60.0f);
    params.push_back(std::move(gainSynthParam));

    auto gainVocParam = std::make_unique<AudioParameterFloat>("gainVoc", "GainVoc", -60.0f, 6.0f, 0.0f);
    params.push_back(std::move(gainVocParam));

    auto lpcVoiceParam = std::make_unique<AudioParameterInt>("lpcVoice", "LpcVoice", 2, 100, 40);
    params.push_back(std::move(lpcVoiceParam));

    auto lpcPitchParam = std::make_unique<AudioParameterInt>("lpcPitch", "LpcPitch", 2, 100, 15);
    params.push_back(std::move(lpcPitchParam));

    auto lpcSynthParam = std::make_unique<AudioParameterInt>("lpcSynth", "LpcSynth", 2, 30, 5);
    params.push_back(std::move(lpcSynthParam));


    auto keyPitchParam =  std::make_unique<AudioParameterChoice>("keyPitch", "KeyPitch",
            StringArray("A", "A#", "B", "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "Chrom"), 12);

    params.push_back(std::move(keyPitchParam));

    params.push_back(std::make_unique<AudioParameterBool>("pitchBool", "PitchBool", true));
    params.push_back(std::make_unique<AudioParameterBool>("vocBool", "VocBool", true));


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
void VocoderAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    ignoreUnused(samplesPerBlock);

    double silenceThresholdDb = -60.0;

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

    // Pitch Corrector
    int corres_256 = floor(256.0 * ratioSR);
    int hopPitch = 3 * corres_256;
    int frameLenPitch = 4 * corres_256;

    pitchProcess.prepare(sampleRate, 100, 800, frameLenPitch, hopPitch, samplesPerBlock, silenceThresholdDb);
    vocoderProcess.prepare(wlenVoc, hopVoc, window_str, silenceThresholdDb);

    int latency = std::max(pitchProcess.getLatency(samplesPerBlock), vocoderProcess.getLatency(samplesPerBlock));
    int samplesToKeep = frameLenPitch;

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

    if (treeState.getRawParameterValue("vocBool")->load())
        vocoderProcess.process(myBuffer);


    if (treeState.getRawParameterValue("pitchBool")->load())
        pitchProcess.process(myBuffer);
    else
        pitchProcess.silence();
    
    auto gainVoice = treeState.getRawParameterValue("gainVoice");
    auto gainSynth = treeState.getRawParameterValue("gainSynth");

    if (gainVoice->load() > -59.0)
        myBuffer.addDryVoice(Decibels::decibelsToGain(gainVoice->load(), -59.0f));

    if (gainSynth->load() > -59.0)
        myBuffer.addSynth(Decibels::decibelsToGain(gainSynth->load(), -59.0f));

    myBuffer.fillOutputBuffer(buffer, nOutputChannels);

}

//==============================================================================
bool VocoderAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

AudioProcessorEditor* VocoderAudioProcessor::createEditor()
{
    return new foleys::MagicPluginEditor (magicState, BinaryData::vocodergui_final, BinaryData::vocodergui_finalSize);
    //return new VocoderAudioProcessorEditor (*this);
}

//==============================================================================
void VocoderAudioProcessor::getStateInformation (MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.

    // MAGIC GUI: let the magicState conveniently handle save and restore the state.
    //            You don't need to use that, but it also takes care of restoring the last editor size
    magicState.getStateInformation (destData);
}

void VocoderAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.

    // MAGIC GUI: let the magicState conveniently handle save and restore the state.
    //            You don't need to use that, but it also takes care of restoring the last editor size
    magicState.setStateInformation (data, sizeInBytes, getActiveEditor());
}

//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new VocoderAudioProcessor();
}
