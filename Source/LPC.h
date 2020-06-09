/*
  ==============================================================================

    LPC.h
    Created: 6 May 2020 11:49:54am
    Author:  Damien Ronssin

  ==============================================================================
*/

#pragma once
#include "MyBuffer.h"
#include "../JuceLibraryCode/JuceHeader.h"
#include <string>
#include <vector>
#include <iostream>
#include <cassert>
#include <cmath>

void lpc(MyBuffer& myBuffer, const double* inputBuffer/*double (MyBuffer::*getSample)(int, int) const*/,
        std::vector<double>& r,
         std::vector<double>& a, std::vector<double>& a_prev, const int& order, const int& wlen,
         const int& startSample, const std::vector<double>& anWindow);

void biaisedAutoCorr(MyBuffer& myBuffer, const double* inputBuffer /*double (MyBuffer::*getSample)(int, int) const*/,
                     std::vector<double>& r, const int& order,  const int& wlen, const int& startSample, const
                     std::vector<double>& anWindow);

void levinsonDurbin(const std::vector<double>& r, std::vector<double>& a, std::vector<double>& a_prev,
                    const int& order);