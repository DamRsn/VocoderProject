/*
  ==============================================================================

    PitchCorrector.cpp
    Created: 6 May 2020 10:24:11am
    Author:  Damien Ronssin

  ==============================================================================
*/

#include "PitchProcess.h"
#include "PluginProcessor.h"


PitchProcess::PitchProcess() {}
PitchProcess::~PitchProcess() {}


int PitchProcess::setLatency() {
    return 0;
}


void PitchProcess::setAudioProcPtr(VocoderAudioProcessor* audioProcPtr)
{
    this->audioProcPtr = audioProcPtr;
}


void PitchProcess::prepare(double fS, double fMin, double fMax, int wlen, int hop, int order, int orderMax,
        double speed)
{
    this->fMin = fMin;
    this->fMax = fMax;
    this->fS = fS;
    this->wlen = wlen;
    this->hop = hop;
    this->order = order;
    this->orderMax = orderMax;
    this->speed = speed;

    this->delta=0.94;
    this->pitch=0;
    this->prevPitch=0;
    this->yinTol = 0.25;

    tauMax = (int)floor(this->fS/this->fMin);
    yinTemp.resize(tauMax, 0.0);


    // Give capacity to vector of pitch mark idx
    int marksPerWindowMax = ceil(wlen * fMax / fS);

    anMarks.reserve(marksPerWindowMax);
    anMarksPrev.reserve(marksPerWindowMax);

    stMarks.reserve(marksPerWindowMax);
    stMarksPrev.reserve(marksPerWindowMax);

    a.resize(orderMax+1, 0.0);
    aPrev.resize(orderMax+1, 0.0);
    r.resize(orderMax+1, 0.0);

    anWindow.resize(wlen, 0.0);
    stWindow.resize(wlen, 1.0);


}


void PitchProcess::process() {

}

int PitchProcess::yin(MyBuffer& myBuffer)
{
    // Returns tau (lag value corresponding to pitch period in number of sample)
    // Update pitch variable
    prevPitch = pitch;
    pitch = 0;

    for (int k; k < tauMax; k++)
    {
        yinTemp[k] = 0;
        for (int i; i < wlen; i++) {
            yinTemp[k] += pow(myBuffer.getVoiceSample(0, i - tauMax)
                                - myBuffer.getVoiceSample(0, i + k - tauMax), 2);
        }
    }

    yinTemp[0] = 1.0;
    double tmp = 0;

    for (int k = 1; k < tauMax; k++)
    {
        tmp += yinTemp[k];
        yinTemp[k] *= k/tmp;
    }

    int tau = floor(fS/fMax);
    while(tau < tauMax)
    {
        if (yinTemp[tau] < yinTol)
        {
            while (yinTemp[tau + 1] < yinTemp[tau])
            {
                tau += 1;
                if (tau+1>=tauMax)
                    break;
            }
            pitch = fS/tau;
            break;
        }
        else
            tau += 1;
    }
    return tau;
}

void PitchProcess::pitchMarks(MyBuffer& myBuffer)
{
    // Copy content of anMarks to anMarksPrev and clear anMarks, there should be no memory allocation
    anMarksPrev = anMarks;
    anMarks.clear();



}

void PitchProcess::getStMarks() {

}


int argExt(MyBuffer& myBuffer, int idxStart, int idxEnd, bool min)
{
    double ext = myBuffer.getVoiceSample(0, idxStart);
    int argExt = 0;

    if (min)
    {
        for (int i = idxStart+1; i < idxEnd; i++)
        {
            if (myBuffer.getVoiceSample(0, i) < ext)
            {
                ext = myBuffer.getVoiceSample(0, i);
                argExt = i;
            }

        }
    }
    else
    {
        for (int i = idxStart+1; i < idxEnd; i++)
        {
            if (myBuffer.getVoiceSample(0, i) > ext)
            {
                ext = myBuffer.getVoiceSample(0, i);
                argExt = i;
            }
        }
    }

    return argExt;
}

