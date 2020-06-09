/*
  ==============================================================================

    Notes.cpp
    Created: 6 May 2020 11:31:57am
    Author:  Damien Ronssin

  ==============================================================================
*/

#include "Notes.h"


Notes::Notes() {}


Notes::~Notes() {}

/**
 * Function to call in prepare to play of pitch process
 * @param key
 * @param fMin : min frequency in Hz
 * @param fMax : max frequency in Hz
 */
void Notes::prepare(Notes::key key, double fMin, double fMax)
{
    freq.reserve(88);
    intervals = {2, 2, 1, 2, 2, 2, 1};

    // Reference frequency
    currKey = key;

    this->fMin = fMin;
    this->fMax = fMax;

    buildFreqVect();
}


/**
 * Build the frequency vector with frequency of given key
 */
void Notes::buildFreqVect()
{
    freq.resize(0);

    // start frequency 30 Hz
    double f = 27.5;
    f = f * pow(2, double(currKey)/12.0);

    int i=0;

    // To avoid compute 2**(1/12) every time
    double factorSemiTone = pow(2, 1.0/12);

    while (freq.empty() || freq.back() < fMax)
    {
        if (currKey != Chrom)
            f = f * pow(factorSemiTone, intervals[i%intervals.size()]);
        else
            f = f * factorSemiTone;

        if (f > fMin)
            freq.push_back(f);

        i+= 1;
    }

    freq.pop_back();
}


/**
 * Return the frequency of the closest note in the given key
 * @param pitch : current pitch in Hz
 * @param key
 * @return
 */
double Notes::getClosestFreq(const double &pitch, Notes::key key) {
    double closestFreq;

    // Check if key has changed, if yes, rebuild freq vect without any memory allocation and update currKey
    if (key != currKey)
    {
        freq.resize(0);
        currKey = key;
        buildFreqVect();
    }

    // use std::lower_bound: gives iterator to first element larger than pitch
    std::vector<double>::iterator lb;
    lb = std::lower_bound(freq.begin(), freq.end(), pitch);

    // get idx of iterator
    int idx = lb - freq.begin();

    if (idx > 0) {
        // Then determine min between idx and idx - 1
        if (abs(freq[idx] - pitch) <= abs(freq[idx - 1] - pitch))
            closestFreq = freq[idx];
        else
            closestFreq = freq[idx - 1];
    }
    else
    {
        closestFreq = freq[idx];
    }

    return closestFreq;
}


