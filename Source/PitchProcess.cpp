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


int PitchProcess::getLatency(int samplesPerBlock) {
    //TODO: something smart that minimizes latency
    int latency;
    if (samplesPerBlock <= frameLen)
        latency = frameLen;

    else
        latency = frameLen;

    return latency;

}


void PitchProcess::setAudioProcPtr(VocoderAudioProcessor* audioProcPtr)
{
    this->audioProcPtr = audioProcPtr;
}


void PitchProcess::prepare(double fS, double fMin, double fMax, int frameLen, int hop, int order, int orderMax,
        double speed, int samplesPerBlock, Notes::key key)
{
    this->fMin = fMin;
    this->fMax = fMax;
    this->fS = fS;
    this->frameLen = frameLen;
    this->hop = hop;
    this->order = order;
    this->orderMax = orderMax;
    this->speed = speed;
    this->samplesPerBlock = samplesPerBlock;

    this->overlap = ((double)(frameLen-hop))/((double)(frameLen));
    //std::cout<< "overlap " << overlap << std::endl;

    this->delta=0.94;
    this->pitch=0;
    this->prevPitch=0;
    this->period=0;
    this->periodNew = 0;
    this->prevVoicedPeriod = 0;
    this->prevPeriod = 0;
    this->beta=1;
    this->yinTol = 0.25;
    this->startSample = 0;
    this->bufferIdxMax = 0;

    this->nChunk = 0;
    this->chunkSize = frameLen - hop;
    this->chunksPerFrame = this->frameLen/chunkSize;


    this->key = key;
    this->outFrame.resize(frameLen, 0.0);

    notes.prepare(key, fMin, fMax);


    tauMax = (int)ceil(this->fS/this->fMin);
    yinTemp.resize(tauMax, 0.0);

    // Give capacity to vector of pitch mark idx
    int marksPerWindowMax = ceil(frameLen * fMax / fS) + 1;

    anMarks.reserve(marksPerWindowMax);
    prevAnMarks.reserve(marksPerWindowMax);
    stMarks.reserve(marksPerWindowMax);
    prevStMarks.reserve(marksPerWindowMax);
    anMarks.resize(0);
    prevAnMarks.resize(0);
    stMarks.resize(0);
    prevStMarks.resize(0);

    a.resize(orderMax+1, 0.0);
    aPrev.resize(orderMax+1, 0.0);
    r.resize(orderMax+1, 0.0);
    // Reserve more memory if future of past samples needed in e
    e.reserve(2*frameLen);
    e.resize(frameLen, 0.0);

    psolaWindow.reserve(2*tauMax + 3);
    periodSamples.reserve(2*tauMax + 3);
    valley = true;

    buildWindows(anWindow, stWindow);
}


void PitchProcess::process(MyBuffer &myBuffer)
{
    if (bufferIdxMax==0)
        bufferIdxMax = myBuffer.getIdxMax();

    while(startSample < myBuffer.getSamplesPerBlock())
    {
        if (nChunk%chunksPerFrame == chunksPerFrame - 1){
            processChunkCont(myBuffer);
            nChunk = 0;
            processChunkStart(myBuffer);
            nChunk += 1;
            nChunk %= chunksPerFrame;
        }

        else if (nChunk == 0) {
            processChunkStart(myBuffer);
            nChunk += 1;
        }
        else {
            processChunkCont(myBuffer);
            nChunk += 1;
        }

        startSample += chunkSize;
    }
    startSample -= myBuffer.getSamplesPerBlock();
}


void PitchProcess::processChunkCont(MyBuffer& myBuffer)
{
    if (!anMarks.empty()){
        psola(myBuffer);
    }

    fillOutputBuffer(myBuffer);
}


void PitchProcess::processChunkStart(MyBuffer& myBuffer)
{
    std::fill(e.begin(), e.end(), 0.0);
    std::fill(outFrame.begin(), outFrame.end(), 0.0);

    yin(myBuffer);
    pitchMarks(myBuffer);
    placeStMarks();

    if (!anMarks.empty()){
        stMarkIdx = 0;
        psola(myBuffer);
    }

    fillOutputBuffer(myBuffer);
}


void PitchProcess::fillOutputBuffer(MyBuffer& myBuffer)
{
    for (int channel = 0; channel < myBuffer.getNumChannels(); channel++) {
        for (int i = 0; i < chunkSize; i++) {
            myBuffer.addOutSample(channel, startSample + i, outFrame[i + nChunk * chunkSize]
                    * stWindow[i + nChunk * chunkSize] * Decibels::decibelsToGain(audioProcPtr->gainPitch));
        }
    }
}


void PitchProcess::yin(MyBuffer& myBuffer)
{
    // Update tau (lag value corresponding to pitch period in number of sample), if no pitch found, tau = tauMax
    // Update pitch variable
    prevPeriod = period;
    prevPitch = pitch;

    if (pitch > 1){
        prevVoicedPeriod = period;
        prevVoicedPitch = pitch;
    }

    pitch = 0;
    period = 0;

    for (int k = 0; k < tauMax; k++){
        yinTemp[k] = 0;
        for (int i = 0; i < frameLen; i++) {
            yinTemp[k] += pow(myBuffer.getVoiceSample(0, startSample + i - tauMax)
                                - myBuffer.getVoiceSample(0, startSample + i + k - tauMax), 2);
        }
    }

    yinTemp[0] = 1.0;
    double tmp = 0;

    for (int k = 1; k < tauMax; k++){
        tmp += yinTemp[k];
        yinTemp[k] *= k/tmp;
    }

    int tau = floor(fS/fMax);

    while(tau < tauMax){
        if (yinTemp[tau] < yinTol){
            while (yinTemp[tau + 1] < yinTemp[tau]){
                tau += 1;
                if (tau+1>=tauMax)
                    break;
            }
            pitch = fS/tau;
            period = tau;
            break;
        }
        else
            tau += 1;
    }
}


void PitchProcess::pitchMarks(MyBuffer& myBuffer)
{
    // Copy content of anMarks to prevAnMarks and clear anMarks, there should be no memory allocation
    prevAnMarks = anMarks;
    anMarks.clear();

    // Number of pitch marks in overlap
    nAnMarksOv = 0;
    prevAnMarks -= hop;

    for (int i=0; i<prevAnMarks.size(); i++){
        if (prevAnMarks[i] >= 0)
            nAnMarksOv += 1;
    }

    bool searchLeft = false;

    // Some ints needed
    int t, l_lim, r_lim, lastMark, sw_c, sw_f;

    // TODO: Check if prevAnMarks empty: use this as a case when there is nothing to do (certain number of unvoiced win)
    // If current frame is voiced
    if (pitch > 1){
        sw_c = floor(delta * period);
        sw_f = ceil((2.0 - delta) * period);

        // If previous frame was also voiced
        if (prevPitch > 1){
            if (nAnMarksOv == 0){
                lastMark = prevAnMarks.back();
                l_lim = std::max(lastMark + std::min(sw_c, int(floor(delta * std::min(prevPeriod, period)))), 0);
                r_lim = std::min(lastMark + std::max(sw_f, int(ceil( (2 - delta)*std::max(prevPeriod, period)))), frameLen);
                t = argExt(myBuffer, l_lim, r_lim, valley);
            } else{
                // Take first of overlapping previous marks
                t = prevAnMarks[prevAnMarks.size() - nAnMarksOv];
            }
        } else {
            // If no info before (previous window was unvoiced): discard any info about previous pitch mark
            // Find arg_ext on all window, and then search right and left
            searchLeft = true;
            t = argExt(myBuffer, 0, frameLen, valley);
        }

        anMarks.push_back(t);

        // Search to the right
        while (anMarks.back() + sw_c < frameLen){
            // Check if all search area is contained in frame
            if (anMarks.back() + sw_f < frameLen)
                anMarks.push_back(argExt(myBuffer, anMarks.back() + sw_c, anMarks.back() + sw_f, valley));
            else{
                if (anMarks.back() + period < frameLen){
                    anMarks.push_back(argExt(myBuffer, anMarks.back() + sw_c, frameLen, valley));
                    break;
                }
                else
                    break;
            }
        }

        // Search to the left if needed
        if (searchLeft){
            while (anMarks.front() - sw_c > 0){
                // Check if all search area is contained in frame
                if (anMarks.front() - sw_f >= 0){
                    anMarks.insert(anMarks.begin(), argExt(myBuffer, anMarks.front() - sw_f,
                            anMarks.front() - sw_c, valley));
                } else {
                    if (anMarks.front() - period >= 0){
                        anMarks.insert(anMarks.begin(), argExt(myBuffer, 0,
                                anMarks.front() - sw_c, valley));
                        break;
                    }
                    else
                        break;
                }
            }
        }
    }
    
    // if current frame is not voiced
    else{
        // if prev marks is not empty (otherwise do nothing so anMarks stays empty)
        if (!prevAnMarks.empty()){
            if (nAnMarksOv > 0){
                for (int i = 0; i < nAnMarksOv; i++)
                    anMarks.push_back(prevAnMarks[prevAnMarks.size() - nAnMarksOv + i]);
            }
            else
                anMarks.push_back(prevAnMarks.back() + prevVoicedPeriod);

            if (prevVoicedPeriod<=0){
                std::cout << "pitchMark: prevVoicedPeriod <= 0" << std::endl;
                assert(false);
            }

            while (anMarks.back() + prevVoicedPeriod < frameLen){
                anMarks.push_back(anMarks.back() + prevVoicedPeriod);
            }
        }
    }
}


void PitchProcess::placeStMarks()
{
    //Todo Check consitence with other update of pitch mark when changing buffer or frame
    prevStMarks = stMarks;
    stMarks.clear();
    nStMarksOv = 0;
    prevStMarks -= hop;

    // Do nothing if no anMarks
    if (anMarks.empty())
        return;

    for (int & prevStMark : prevStMarks){
        if (prevStMark >= 0)
            nStMarksOv += 1;
    }

    int firstMark;
    prevClosestFreq = closestFreq;

    if (pitch > 1){
        closestFreq = notes.getClosestFreq(pitch, key);
        beta = closestFreq / pitch;
        periodNew = round(period/beta);
    }
    else{
        closestFreq = 0;
        periodNew = prevVoicedPeriod;
    }

    if (periodNew <= 0)
    {
        std::cout << "placeStMarks: periodNew0 <= 0" << std::endl;
        assert(false);
    }

    // Todo: deal with if the closest note has changed or not
    // Todo: change this mess
    // If curr frame voiced
    if (pitch > 1){
        // if previous frame also voiced
        if (prevPitch > 1){
            if (nStMarksOv > 0)
                firstMark = prevStMarks[prevStMarks.size() - nStMarksOv];
            else
                if (prevStMarks.back() + periodNew>=0)
                    firstMark = prevStMarks.back() + periodNew;
                else
                    firstMark = anMarks[0];
        }
        // if previous frame was unvoiced
        else
            firstMark = anMarks[0];
    }
    // if current frame unvoiced
    else{
        if (prevStMarks.empty())
            return;

        if (nStMarksOv > 0)
            firstMark = prevStMarks[prevStMarks.size() - nStMarksOv];
        else{
            if (periodNew<=0){
                std::cout << "placeStMarks: periodNew <= 0" << std::endl;
                assert(false);
            }

            int n = 1;
            while(prevStMarks.back() + n*periodNew < 0)
                n +=1;

            firstMark = prevStMarks.back() + n * periodNew;
        }
    }

    stMarks.push_back(firstMark);

    if (periodNew <= 0){
        std::cout << "placeStMarks: periodNew2 <= 0" << std::endl;
        assert(false);
    }

    while(stMarks.back() + periodNew < frameLen)
        stMarks.push_back(stMarks.back() + periodNew);

}


void PitchProcess::psola(MyBuffer& myBuffer)
{
    int periodPsola;

    if (pitch > 1)
        periodPsola = period;
    else
        periodPsola = prevVoicedPeriod;

    // cl stands for closest
    int stMark, clIdx, clAnMark, startIdx, stopIdx;
    periodSamples.resize(2*periodPsola+1);
    xInterp.resize(2*periodPsola+1);
    fillPsolaWindow(psolaWindow, period);

    while(stMarkIdx < stMarks.size()){
        stMark = stMarks[stMarkIdx];

        // check if we have processed all the necessary stMarks for this buffer
        // TODO: Here to decide where to stop for a chunk
        if (stMark - periodPsola >= (nChunk + 1) * chunkSize)
            break;

        // Return closest complete anMarks (can be in previous)
        clIdx = getClosestAnMarkIdx(anMarks, prevAnMarks, stMark, periodPsola);

        if (clIdx >= 0)
            clAnMark = anMarks[clIdx];
        else
            clAnMark = prevAnMarks[prevAnMarks.size() - clIdx];

        if (stMarkIdx !=0 && stMarkIdx != stMarks.size()-1){
            for (int j = 0; j < 2 * periodPsola + 1; j++) {
                periodSamples[j] = myBuffer.getVoiceSample(0,
                        - nChunk * chunkSize + startSample + clAnMark - periodPsola + j)
                        * psolaWindow[j];
                xInterp[j] = stMark + (-periodPsola + j) / beta;
            }

            startIdx = std::max(int(floor(xInterp[0])), 0);
            stopIdx = std::min(int(ceil(xInterp.back())), frameLen);

            interp(xInterp, periodSamples, outFrame, psolaWindow, startIdx, stopIdx);
        }else if (stMarkIdx == 0){
            for (int j = 0; j < 2 * periodPsola + 1; j++) {
                if (j < periodPsola)
                    periodSamples[j] = myBuffer.getVoiceSample(0,
                                       - nChunk * chunkSize  + startSample + clAnMark - periodPsola + j);
                else
                    periodSamples[j] = myBuffer.getVoiceSample(0,
                                        - nChunk * chunkSize + startSample + clAnMark - periodPsola + j)
                            * psolaWindow[j];

                xInterp[j] = stMark + (- periodPsola + j) / beta;
            }

            startIdx = std::max(int(floor(xInterp[0])), 0);
            stopIdx =  std::min(int(ceil(xInterp.back())), frameLen);

            interp(xInterp, periodSamples, outFrame, psolaWindow, startIdx, stopIdx);
        }else{
            for (int j = 0; j < 2 * periodPsola + 1; j++) {
                if (j < periodPsola)
                    periodSamples[j] = myBuffer.getVoiceSample(0,
                               - nChunk * chunkSize  + startSample + clAnMark - periodPsola + j)
                            * psolaWindow[j];
                else
                    periodSamples[j] = myBuffer.getVoiceSample(0,
                            - nChunk * chunkSize + startSample + clAnMark - periodPsola + j);

                xInterp[j] = stMark + (- periodPsola + j) / beta;
            }

            startIdx = std::max(int(floor(xInterp[0])), 0);
            stopIdx = std::min(int(ceil(xInterp.back())), frameLen);

            interp(xInterp, periodSamples, outFrame, psolaWindow, startIdx, stopIdx);
        }

        stMarkIdx += 1;
    }
}


int PitchProcess::argExt(const MyBuffer& myBuffer, int idxStart, int idxEnd, bool min)
{
    double ext = myBuffer.getVoiceSample(0, startSample + idxStart);
    int argExt = idxStart;

    if (min){
        for (int i = idxStart+1; i < idxEnd; i++){
            if (myBuffer.getVoiceSample(0, startSample + i) < ext){
                ext = myBuffer.getVoiceSample(0, startSample + i);
                argExt = i;
            }
        }
    }
    else
    {
        for (int i = idxStart+1; i < idxEnd; i++){
            if (myBuffer.getVoiceSample(0, startSample + i) > ext){
                ext = myBuffer.getVoiceSample(0, startSample + i);
                argExt = i;
            }
        }
    }

    return argExt;
}


int PitchProcess::getClosestAnMarkIdx(const std::vector<int>& anMarks, const std::vector<int>& prevAnMarks,
        const int& stMark, int periodPsola)
{
    int closestIdx;
    auto lb = std::lower_bound(anMarks.begin(), anMarks.end(), stMark);

    // get idx of iterator
    int idx = lb - anMarks.begin();

    if (idx > 0 && idx < anMarks.size()) {
        // Then determine min between idx and idx - 1 and that the chosen one is complete
        if (abs(anMarks[idx] - stMark) <= abs(anMarks[idx - 1] - stMark) &&
            anMarks[idx] + periodPsola - nChunk * chunkSize < bufferIdxMax - startSample)
            closestIdx = idx;

        else if (anMarks[idx - 1] + periodPsola - nChunk * chunkSize < bufferIdxMax - startSample)
            closestIdx = idx - 1;

        else
            if (idx - 2 > 0)
                closestIdx = idx - 2;

            // Take last element of non overlapping of prevAnMarks
            else
                closestIdx = - nAnMarksOv - 1;
    }
    else if (idx==0)
        closestIdx = idx;

    else if (idx==anMarks.size()){
        if (anMarks[idx] + periodPsola - nChunk * chunkSize < bufferIdxMax - startSample)
            closestIdx = idx - 1;
        else
            if (idx - 2 >=0)
                closestIdx = idx-2;
            else
                assert(false);
    }

    else
    {
        std::cout << "anMarks= " << anMarks;
        std::cout << "idx = " << idx << std::endl;
        std::cout << "Error in getClosestAnMarkIdx"<< std::endl;
        assert(false);
    }

    return closestIdx;
}


void PitchProcess::interp(vector<double> &x, const vector<double> &y, std::vector<double>& outFrame,
        const vector<double>& psolaWindow, const int& startIdx, const int &stopIdx)
{
    auto startSearchIt = x.begin();
    std::vector<double>::iterator lb;
    int lbIdx;
    double value;
    //std::cout << "in interp";
    for (int i = startIdx; i < stopIdx; i++){
        if (i >= x.front() && i <= x.back()){
            // Find interval
            lb = std::lower_bound(startSearchIt, x.end(), i);

            // update startsearch iterator (because we know that everything is ordered)
            startSearchIt = lb - 1;
            lbIdx = lb - x.begin();

            if (lbIdx>0)
                value = y[lbIdx - 1] + (y[lbIdx] - y[lbIdx-1])/(x[lbIdx] - x[lbIdx-1]) * (i - x[lbIdx-1]);
            else
                value = y[lbIdx];

            outFrame[i] += value;// * psolaWindow[i - startIdx];
        }

        else if (i > x.back())
            break;
    }
}


void PitchProcess::fillPsolaWindow(std::vector<double>& psolaWindow, const int& T) {
    psolaWindow.resize(2*T+1);
    dsp::WindowingFunction<double>::fillWindowingTables(&psolaWindow[0], 2*T+1,
                                                        dsp::WindowingFunction<double>::hann,false);
}


void PitchProcess::buildWindows(std::vector<double>& anWindow, std::vector<double>& stWindow)
{
    anWindow.resize(frameLen, 1.0);
    stWindow.reserve(frameLen);
    stWindow.resize(2*int(round(overlap*frameLen)));

    dsp::WindowingFunction<double>::fillWindowingTables(&stWindow[0], stWindow.size(),
                                                        dsp::WindowingFunction<double>::hann,false);

    for (int i=0; i < int(frameLen - 2*int(round(overlap*frameLen))); i++)
        stWindow.insert(stWindow.begin() + int(round(overlap*frameLen)), 1.0);

    if (stWindow.size() != frameLen){
        std::cout << "stWindow.size() != frameLen" << std::endl;
        assert(false);
    }
}



bool equal(double a, double b)
{
    if (abs(a-b)<pow(10, -12))
        return true;
    else
        return false;
}


template <typename T>
void operator-=(std::vector<T>& v1, const T& a)
{
    for(int & el : v1)
        el -= a;
}

template <typename T>
void operator+=(std::vector<T>& v1, const T& a)
{
    for(int & el : v1)
        el += a;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
    os << "[";
    for (int i = 0; i < v.size(); ++i) {
        os << v[i];
        if (i != v.size() - 1)
            os << ", ";
    }
    os << "]\n";
    return os;
}











