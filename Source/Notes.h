/*
  ==============================================================================

    Notes.h
    Created: 6 May 2020 11:31:57am
    Author:  Damien Ronssin

  ==============================================================================
*/

#pragma once
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iostream>



class Notes
{
public:
    Notes();
    ~Notes();

    enum key {A=0, As=1, Bb=1, B=2, C=3, Cs=4, Db=4, D=5, Ds=6, Eb=6, E=7, F=8, Fs=9, Gb=9, G=10, Gs=11, Ab = 11, Chrom=12};
    void prepare(key key, double fMin, double fMax);
    double getClosestFreq(const double& pitch, key key);

private:
    double fMin;
    double fMax;

    // Array with all frequencies from 60 hz to 1920 Hz (1920 not in)
    std::vector<double> freq;

    // Intervals of a major key
    std::vector<int> intervals;
    void buildFreqVect();

    key currKey;
};
