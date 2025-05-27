Following scripts are custom made for Matlab and are made to extract data from patch clamp experiment traces acquired with Igor.
Patch clamp experiment acquire with Igor give .ibw files.
With the following script we use the function IBWRead to convert .ibw files to .mat files 

Each trace is a separate file within its cell folder. 

To use the scripts:

First save all the given function in a folder and add it to your matlab path
The main script is called IgorToMatlab_AnalysesGH.m, all the functions are implemented within it. 

The script extract each trace of a patch cell and extract the following parameters when available:

- Firing frequency
- Capacitance
- Input resistance
- Tau
- Resting membrane potential
- Rheobase
- Rise Time, duration and speed
- Amplitude of action potential
- Decay Time, duration and speed
- Half width and half width duration
- AHP Time, amplitude and speed
- DAP
- 1st AP delay
- interspike interval
- adapataion index
- accomodation index
- SAG
