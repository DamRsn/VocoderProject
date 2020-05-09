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


int PitchProcess::getLatency() {
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
    int marksPerWindowMax = ceil(wlen * fMax / fS) + 1;

    anMarks.reserve(marksPerWindowMax);
    prevAnMarks.reserve(marksPerWindowMax);

    stMarks.reserve(marksPerWindowMax);
    prevStMarks.reserve(marksPerWindowMax);

    a.resize(orderMax+1, 0.0);
    aPrev.resize(orderMax+1, 0.0);
    r.resize(orderMax+1, 0.0);

    anWindow.resize(wlen, 0.0);
    stWindow.resize(wlen, 1.0);

    psolaWindow.reserve(2*tauMax + 3);
    valley = true;
}


void PitchProcess::process() {

}

int PitchProcess::yin(MyBuffer& myBuffer)
{
    // Returns tau (lag value corresponding to pitch period in number of sample), if no pitch found, tau = tauMax
    // Update pitch variable

    prevPeriod = period;
    prevPitch = pitch;
    pitch = 0;

    for (int k = 0; k < tauMax; k++)
    {
        yinTemp[k] = 0;
        for (int i = 0; i < wlen; i++) {
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

    // update prevVoidPeriod
    if (pitch > 1)
        prevVoicedPeriod = tau;

    return tau;
}


void PitchProcess::pitchMarks(MyBuffer& myBuffer)
{
    // Copy content of anMarks to prevAnMarks and clear anMarks, there should be no memory allocation
    prevAnMarks = anMarks;
    anMarks.clear();

    // Number of pitch marks in overlap
    int nMarksOv = 0;
    for (int i=0; i<prevAnMarks.size(); i++)
    {
        prevAnMarks[i] -= hop;
        if (prevAnMarks[i] >= 0)
            nMarksOv += 1;
    }

    bool searchLeft = false;

    // Some ints needed
    int t, l_lim, r_lim, lastMark, sw_c, sw_f;

    // TODO: Check if prevAnMarks empty: use this as a case when there is nothing to do (certain number of unvoiced win)
    // If current frame is voiced
    if (pitch > 1)
    {
        sw_c = floor(delta * period);
        sw_f = ceil(delta * period);

        // If previous frame was also voiced
        if (prevPitch > 1)
        {
            if (nMarksOv == 0)
            {
                lastMark = prevAnMarks.back();
                l_lim = std::max(lastMark + std::min(sw_c, int(floor(delta * std::min(prevPeriod, period)))), 0);
                r_lim = std::min(lastMark + std::max(sw_f, int(ceil( (2 - delta)*std::max(prevPeriod, period)))), wlen);
                t = argExt(myBuffer, l_lim, r_lim, valley);
            } else{
                // Take first of overlapping previous marks
                t = prevAnMarks[prevAnMarks.size() - nMarksOv];
            }
        } else {
            // If no info before (previous window was unvoiced): discard any info about previous pitch mark
            // Find arg_ext on all window, and then search right and left
            searchLeft = true;
            t = argExt(myBuffer, 0, wlen, valley);
        }

        anMarks.push_back(t);

        // Search to the right
        while (anMarks.back() + sw_c < wlen)
        {
            // Check if all search area is contained in frame
            if (anMarks.back() + sw_f < wlen)
                anMarks.push_back(argExt(myBuffer, anMarks.back() + sw_c, anMarks.back() + sw_f, valley));
            else
            {
                if (anMarks.back() + period < wlen)
                {
                    anMarks.push_back(argExt(myBuffer, anMarks.back() + sw_c, wlen, valley));
                    break;
                }
                else
                    break;
            }
        }

        // Search to the left if needed
        if (searchLeft)
        {
            while (anMarks.back() - sw_c > 0)
            {
                // Check if all search area is contained in frame
                if (anMarks.back() - sw_f >= 0)
                {
                    anMarks.insert(anMarks.begin(), argExt(myBuffer, anMarks.back() - sw_f,
                            anMarks.back() - sw_c, valley));
                } else {
                    if (anMarks.back() - period >= 0)
                    {
                        anMarks.insert(anMarks.begin(), argExt(myBuffer, 0,
                                anMarks.back() - sw_c, valley));
                        break;
                    }
                    else
                        break;
                }
            }
        }
    }
    
    // if current frame is not voiced
    else
    {
        // if prev marks is not empty (otherwise do nothing so anMarks stays empty)
        if (!prevAnMarks.empty())
        {
            if (nMarksOv > 0)
            {
                for (int i = 0; i < nMarksOv; i++)
                    anMarks.push_back(prevAnMarks[prevAnMarks.size() - nMarksOv + i]);

            }
            else
                anMarks.push_back(prevAnMarks.back() + prevVoicedPeriod);
            
            while (anMarks.back() + prevVoicedPeriod)
                anMarks.push_back(anMarks.back() + prevVoicedPeriod);
        }
    }


}

void PitchProcess::placeStMarks()
{
    prevStMarks = stMarks;
    stMarks.clear();
    int nMarksOv = 0;

    for (int & prevAnMark : prevAnMarks)
    {
        prevAnMark -= hop;
        if (prevAnMark >= 0)
            nMarksOv += 1;
    }

    int firstMark;
    prevClosestFreq = closestFreq;

    if (pitch > 1)
    {
        closestFreq = notes.getClosestFreq(pitch, key);
        beta = closestFreq / pitch;
        periodNew = round(period/beta);
    }
    else
    {
        closestFreq = 0;
        beta = 0;
        periodNew = prevVoicedPeriod;
    }


    // Todo: deal with if the closest note has changed or not

    // Todo change this mess
    // If curr frame voiced
    if (pitch > 1)
    {
        // if previous frame also voiced
        if (prevPitch > 1)
        {
            if (nMarksOv > 0)
                firstMark = prevStMarks[prevStMarks.size() - nMarksOv];
            else
                if (prevStMarks.back() + periodNew>=0)
                    firstMark = prevStMarks.back() + periodNew;
                else
                    firstMark = anMarks[0];
        }
        // if previous frame was unvoiced
        else
        {
            firstMark = anMarks[0];
        }
    }
    // if current frame unvoiced
    else
    {
        if (prevStMarks.empty())
            return;

        if (nMarksOv > 0)
            firstMark = prevStMarks[prevStMarks.size() - nMarksOv];
        else
            firstMark = prevStMarks.back() + periodNew;
    }


    stMarks.push_back(firstMark);

    while(stMarks.back() + periodNew < wlen)
        stMarks.push_back(stMarks.back() + periodNew);

}


void PitchProcess::fillPsolaWindow(std::vector<double>& psolaWindow, const int& T) {
    psolaWindow.resize(2*T+1);
    dsp::WindowingFunction<double>::fillWindowingTables(&psolaWindow[0], 2*T+1,
                                                        dsp::WindowingFunction<double>::hann,false);
}


int PitchProcess::argExt(const MyBuffer& myBuffer, int idxStart, int idxEnd, bool min)
{
    double ext = myBuffer.getVoiceSample(0, idxStart);
    int argExt = idxStart;

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

bool equal(double a, double b)
{
    if (abs(a-b)<pow(10, -12))
        return true;
    else
        return false;
}