# STFT Phase Reconstruction for Speech Enhancement

The algorithm estimates the phase spectrum of the underlying speech signal. The voiced part of the speech is modeled using the harmonic model. As a prerequisite step, the fundamental frequencies of the voiced speech are estimated. The phase reconstruction algorithm makes it possible to enhance speech using only the fundamental frequencies and the noisy signal.

I later improved upon this algorithm by using pitch-synchronous window size for STFT instead of constant window size. This made the estimation of phase spectrum faster and more accurate. Visit the [project's repository](https://github.com/geetkhatri/speech-enhancement-psr) to see how it works.

### References

The MATLAB code for the implementation of the algorithm makes use of:

- Mike Brookes. VOICEBOX: A speech processing toolbox for MATLAB (http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html).
- IoSR Matlab Toolbox (https://github.com/IoSR-Surrey/MatlabToolbox).
