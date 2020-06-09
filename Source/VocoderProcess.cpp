    /*
  ==============================================================================

    VocoderProcess.cpp
    Created: 20 Apr 2020 10:17:42am
    Author:  Damien Ronssin

  ==============================================================================
*/

#include "VocoderProcess.h"
#include "PluginProcessor.h"
#define PI 3.14159265


VocoderProcess::VocoderProcess(){}

double sum(const std::vector<double>& v);
double maxAbs(const std::vector<double>& v);

void shift(std::vector<double>& v);

void VocoderProcess::setAudioProcPtr(VocoderAudioProcessor* audioProcPtr)
{
    this->audioProcPtr = audioProcPtr;
}

/**
 * Function to call in prepare to play to allocate memory and instantiate variables for the vocoder
 * @param wlen
 * @param hop
 * @param windowType
 * @param silenceThresholdDb
 */
void VocoderProcess::prepare(int wlen, int hop, std::string windowType, double silenceThresholdDb)
{
    this->wlen = wlen;
    this->hop = hop;
    startSample = 0;

    this->orderVoice = audioProcPtr->treeState.getRawParameterValue("lpcVoice")->load();
    this->orderMaxVoice = audioProcPtr->treeState.getParameterRange("lpcVoice").end;


    this->orderSynth =  audioProcPtr->treeState.getRawParameterValue("lpcSynth")->load();
    this->orderMaxSynth = audioProcPtr->treeState.getParameterRange("lpcSynth").end;
    this->silenceThresholdDb = silenceThresholdDb;

    rVoice.resize(orderMaxVoice + 1, 1.0);
    aVoice.resize(orderMaxVoice + 1, 1.0);
    aPrevVoice.resize(orderMaxVoice + 1, 1.0);
    eVoice.resize(wlen, 0.0);

    rSynth.resize(orderMaxSynth + 1, 1.0);
    aSynth.resize(orderMaxSynth + 1, 1.0);
    aPrevSynth.resize(orderMaxSynth + 1, 1.0);
    eSynth.resize(wlen, 0.0);

    out.resize(wlen, 0.0);

    // gain business
    lambda = 0.0;
    g = 0.0;
    EeSynth = 1.0;
    EeVoice = 0.0;

    EeVoiceArr.resize(10, 0.0);
    EeSynthArr.resize(10, 0.0);

    setWindows(windowType);
}

VocoderProcess::~VocoderProcess() {}

/**
 * return latency needed for the vocoder to work
 * @param samplesPerBlock : buffer size
 * @return
 */
int VocoderProcess::getLatency(int samplesPerBlock)
{
    // TODO: can do better?
    int latency;

    // Bypass everything: set latency to wlen so we are sure no pb since samples per block > 0
    latency = wlen;
    return latency;
}


/**
 * Resize and fill anWindow and stWindow
 * @param windowType : "hann" or "sine"
 */
void VocoderProcess::setWindows(std::string windowType)
{
    // Build window vector for analysis and synthesis

    anWindow.resize(wlen, 0.0);
    stWindow.resize(wlen, 0.0);

    double overlap = double((wlen-hop))/double(wlen);

    double overlapFactor = 1.0;

    if (abs(overlap - 0.75) < pow(10, -10)){
        overlapFactor = 1.0/sqrt(2);
    }

    if (abs(overlap - 0.75) > pow(10, -10) &&
        abs(overlap - 0.5) > pow(10, -10)){
        std::cout << "Invalid overlap" << std::endl;
        assert(false);
    }

    if (windowType=="hann"){
        dsp::WindowingFunction<double>::fillWindowingTables(&stWindow[0], wlen,
                                                           dsp::WindowingFunction<double>::hann,false);
        for (int i = 0; i < wlen; i++)
            stWindow[i] *= overlapFactor;


        std::fill(anWindow.begin(), anWindow.end(), 1.0);
    }
    else if (windowType=="sine"){
        for (int i = 0; i < wlen; i++){
            anWindow[i] = overlapFactor * sin((i+0.5)*PI/double(wlen));
            stWindow[i] = overlapFactor * sin((i+0.5)*PI/double(wlen));
        }
    }
    else{
        std::cerr << "Unknown window type" << std::endl;
        assert(false);
    }
}


/**
 * Change the LPC order for voice. Reads the value from value tree state
 */
void VocoderProcess::setOrderVoice()
{
    auto orderVoicePtr = audioProcPtr->treeState.getRawParameterValue("lpcVoice");
    orderVoice = orderVoicePtr->load();
    if (orderVoice > orderMaxVoice){
        std::cout << "New order for LPC Voice is larger than OrderMax" << std::endl;
        assert(false);
    }

}


/**
 * Change the LPC order for synth. Reads the value from value tree state
 */
void VocoderProcess::setOrderSynth()
{
    auto orderSynthPtr = audioProcPtr->treeState.getRawParameterValue("lpcSynth");
    orderSynth = orderSynthPtr->load();

    if (orderSynth > orderMaxSynth){
        std::cout << "New order for LPC Synth is larger than OrderMax" << std::endl;
        assert(false);
    }
}


/**
 * Function called in processBlock. Launch different functions at appropriate moments. Update the class counter
 * startSample
 * @param myBuffer
 */
void VocoderProcess::process(MyBuffer &myBuffer)
{
    // Call process window
    while (startSample < myBuffer.getSamplesPerBlock()){
        processWindow(myBuffer);
        startSample += hop;
    }

    // Update startSample
    startSample -= myBuffer.getSamplesPerBlock();
}


/**
 * Process a window (frame). Function called by process
 * @param myBuffer
 */
void VocoderProcess::processWindow(MyBuffer &myBuffer)
{
    // Update order of LPC
    setOrderSynth();
    setOrderVoice();

    //double RMSVoiceDb = Decibels::gainToDecibels(myBuffer.getRMSLevelVoice(startSample, wlen));
    //double RMSSynthDb = Decibels::gainToDecibels(myBuffer.getRMSLevelSynth(startSample, wlen));

    double RMSVoiceDb = Decibels::gainToDecibels(myBuffer.getRMSLevelVoiceFull());
    double RMSSynthDb = Decibels::gainToDecibels(myBuffer.getRMSLevelSynthFull());


    if (RMSVoiceDb < silenceThresholdDb or RMSSynthDb < silenceThresholdDb)
        return;

    // LPC on Voice
    lpc(myBuffer, myBuffer.getVoiceReadPtr(), rVoice, aVoice, aPrevVoice, orderVoice,
            wlen, startSample, anWindow);

    // LPC on Synth
    lpc(myBuffer, myBuffer.getSynthReadPtr(0), rSynth, aSynth, aPrevSynth, orderSynth,
            wlen, startSample, anWindow);

    // Get excitation signal on voice, with energy
    filterFIR(myBuffer, &MyBuffer::getVoiceSample, eVoice, aVoice, orderVoice, EeVoice);

    // Get excitation signal on synth, with energy
    filterFIR(myBuffer, &MyBuffer::getSynthSample, eSynth, aSynth, orderSynth, EeSynth);

    // Final: get final cross-synthesis signal: excitation of synth and enveloppe of voice
    filterIIR(myBuffer, aVoice, orderVoice);

}


/**
 * Compute the residual for voice or synth
 * @param myBuffer
 * @param getSample : pointer to function getVoiceSample or getSynthSample
 * @param e : residual vector
 * @param a : LPC filter coeff
 * @param order : LPC order
 * @param E : energy of residual
 */
void VocoderProcess::filterFIR(MyBuffer& myBuffer, double (MyBuffer::*getSample)(int, int) const,
        std::vector<double>& e, const std::vector<double>& a, const int order, double& E)
{
    // Get excitation signal and its energy on current window
    E = 0.0;
    for (int i = 0; i < wlen; i++){
        e[i] = a[0] * (myBuffer.*getSample)(0, startSample + i) * anWindow[i];

        for (int k = 1; k < order + 1; k++){
            if (i - k >= 0)
                e[i] += (myBuffer.*getSample)(0, startSample + i - k) * anWindow[i - k] * a[k];
            else
                break;
        }
    E += e[i]*e[i];
    }
}


/**
 * IIR filter and add output samples to ouput buffer of myBuffer
 * @param myBuffer
 * @param a : LPC filter coeff of voice
 * @param order : order of voice LPC filter
 */
void VocoderProcess::filterIIR(MyBuffer& myBuffer, const std::vector<double>& a,
        const int order)
{
    // Array of energy of previous windows, use mean to have smooth gain changes
    shift(EeVoiceArr);
    shift(EeSynthArr);
    EeVoiceArr[0] = EeVoice;
    EeSynthArr[0] = EeSynth;

    // TODO: this is bad if using something that is not a synth with lots of power, can do better
    if (EeSynth > pow(10, -4))
        // For now lambda is 0
        g = sqrt(sum(EeVoiceArr)/sum(EeSynthArr));

    else
        g = 0.0;

    for (int i = 0; i < wlen; i++){
        out[i] =  g * eSynth[i];

        for (int k = 1; k < order + 1; k++) {
            if (i-k >= 0)
                out[i] -= out[i-k] * a[k];
            else
                break;
        }
    }

    auto gainVoc = audioProcPtr->treeState.getRawParameterValue("gainVoc");

    // Add samples to output buffer
    for (int channel = 0; channel < myBuffer.getNumOutChannels(); channel++){
        for (int i = 0; i < wlen; i++)
            myBuffer.addOutSample(channel, startSample + i,
                    Decibels::decibelsToGain(gainVoc->load(), -59.0f) * out[i] * stWindow[i]);

    }
}


// Some methods on vector
double sum(const std::vector<double>& v)
{
    double s = 0;
    for (auto& n : v)
        s += n;
    return s;
}


double maxAbs(const std::vector<double>& v)
{
    double max = 0.0;
    for (auto& n : v){
        if (abs(n) > max)
            max = abs(n);
    }
    return max;
}


void shift(std::vector<double>& v)
{
    for (int i = v.size()-1; i>0; i--)
        v[i]=v[i-1];

    v[0]=0.0;
}



