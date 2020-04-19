/*
  ==============================================================================

    SynthSound.h
    Created: 15 Apr 2020 11:52:29am
    Author:  Damien Ronssin

  ==============================================================================
*/

#pragma once
#include "../JuceLibraryCode/JuceHeader.h"

class SynthSound: public SynthesiserSound
{
public:
    bool appliesToNote(int /*midinoteNumbers*/)
    {
        return true;
    }
    
    bool appliesToChannel(int /*midiChannel*/)
    {
        return true;
    }
};
