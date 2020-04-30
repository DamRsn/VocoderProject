/*
  ==============================================================================

    SynthVoice.cpp
    Created: 28 Apr 2020 1:28:06pm
    Author:  Damien Ronssin

  ==============================================================================
*/

#include "SynthVoice.h"
#include "PluginProcessor.h"

/*
SynthVoice::SynthVoice(VocoderAudioProcessor* audioProcPtr)
{
    this->audioProcPtr = audioProcPtr;
}
 */

bool SynthVoice::canPlaySound (SynthesiserSound* sound)
{
    return dynamic_cast <SynthSound*>(sound) != nullptr;
}

void SynthVoice::startNote (int midiNoteNumber, float velocity, SynthesiserSound *sound, int currentPitchWheelPosition)
{
    env1.trigger = 1;

    // Todo: deal with velocity in a better way
    //level = audioProcPtr->velocity;
    level = 0.1;
    frequency = MidiMessage::getMidiNoteInHertz(midiNoteNumber);
}

void SynthVoice::stopNote (float velocity, bool allowTailOff)
{
    env1.trigger = 0;
    allowTailOff = true;

    if (velocity == 0)
        clearCurrentNote();
}

void SynthVoice::pitchWheelMoved (int newPitchWheelValue)
{

}

void SynthVoice::controllerMoved (int controllerNumber, int newControllerValue)
{

}

void SynthVoice::renderNextBlock (AudioBuffer<float> &outputBuffer, int startSample, int numSamples)
{
    env1.setAttack(200);
    env1.setDecay(200);
    env1.setSustain(0.9);
    env1.setRelease(100);

    for (int sample = 0; sample < numSamples; sample++)
    {
        double theWave = osc1.saw(frequency);
        double theSound = env1.adsr(theWave, env1.trigger) * level;

        for (int channel = 0; channel < outputBuffer.getNumChannels(); ++channel)
        {
            outputBuffer.addSample(channel, sample, theSound);
        }
        ++startSample;
    }
}
