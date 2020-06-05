# VocoderProject
Juce C++ code and other related material for EPFL Semester Project. Implementation of a Vocoder + Pitch corrector audio plugin.
The plugin is still beta version. Read *Use the plugin* section.
A presentation video is available [here](https://drive.google.com/open?id=1OAemdORaguzs_O-0Av5Q_HthpCF036fL "Presentation video").

### Source
All the juce code is available in the folder Source/. 
The code will be properly commented soon.

### Vocoder.jucer
File for the projucer. 
Note that the module gui magic by ffaudio is needed available [here](https://github.com/ffAudio/foleys_gui_magic "Foley's GUI magic").

### Notebook
Interactive presentation of the algorithms in a jupyter notebook.

### SoundSamples
Some sound samples to listen obtained with the plugin.

### Ressources 
xml file of the user interface

### Plugin
VST and AU files of the plugin. (beta version !)

## Use the plugin
Download the file corresponding to the plugin type your DAW needs and place it in your Plugin folder. 
The folder path depend on your OS. Then import the plugin to your DAW (a plugin scan may be needed) and
add it to an audio track. Like this, the pitch corrector should work. For the vocoder, you should
redirect the sound of your synthesizer to the plugin via Side-chain. The sidechain configuration
process depends on the DAW but it should be easy as a lot of commercial plugins use sidechaining
to get multiple audio inputs. Note that the plugins receives only mono for the main input (voice),
stereo for the side chain input (synthesizer) and stereo for output.
As said before, the current plugin is still a beta version and has not been fully tested in all DAWs
and with all sampling rates. Some problems may still occur. Start using it with low volume and
maximum buffer size to check if it works properly.
Only 2 DAWS were tested (Logic Pro X and REAPER) and only on my own computer. For
REAPER the plugin works for all buffer size above 128 samples and for Logic it works only with
the largest buffer size 1024. The plugin has only been tested with a sampling rate of 44100 Hz. The
plugin has been implemented so that it can work with any sampling rate and produce the same
results, but since no test has been performed, I can't ensure that it works.
