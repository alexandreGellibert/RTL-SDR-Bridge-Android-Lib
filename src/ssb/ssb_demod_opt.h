//
// Created by florent on 13/11/2025.
//

#ifndef SSB_DEMOD_OPT_H
#define SSB_DEMOD_OPT_H

#include <cstdint>
#include <vector>
#include <complex>


// Simple IIR 2nd order helper
struct IIR2 {
    float a0=0, a1=0, a2=0;
    float b1=0, b2=0;
    float z1=0, z2=0;
};


// Standard Biquad
struct Biquad {
    float a0=0, a1=0, a2=0;
    float b1=0, b2=0;
    float z1=0, z2=0;
};

// Complex (I/Q-preserving) lowpass: unlike iir2Process (which filters only the
// real channel and collapses I/Q into a real-only signal — fine for the SSB
// phasing detector but destroys the phase relationship an FM discriminator
// needs), this filters the real and imaginary channels independently so the
// complex IQ signal stays intact.
struct IIR2Complex {
    IIR2 real;
    IIR2 imag;
};

// ----------------------------------------------------
// Fonctions de traitement du pipeline SSB
// ----------------------------------------------------

// 1) Conversion IQ 8-bit -> complex<float>
void convertIQ(const unsigned char* buffer, int len, std::vector<std::complex<float>>& iq);

// 2) Suppression de la composante DC
void removeDC(std::vector<std::complex<float>>& iq, float alpha);

// 3) Filtre passe-bas IIR2
void iir2InitLowpass(IIR2 &f, float fs, float fc, float Q);
void iir2Process(IIR2 &f, std::vector<std::complex<float>> &sig);

// 3b) Complex-preserving lowpass (see IIR2Complex above)
void iir2InitLowpassComplex(IIR2Complex &f, float fs, float fc, float Q);
void iir2ProcessComplex(IIR2Complex &f, std::vector<std::complex<float>> &sig);

// 4) Démodulation SSB
void demodSSB(const std::vector<std::complex<float>>& iq, std::vector<float>& audio, bool upper);

// 4b) FM discriminator: instantaneous phase difference between consecutive
// samples, arg(cur * conj(prev)), normalized to roughly [-1, 1]. Requires the
// real complex IQ signal (use iir2ProcessComplex upstream, not iir2Process).
void demodFM(const std::vector<std::complex<float>>& iq, std::vector<float>& audio);

// 5) AGC amélioré
void adaptiveAGC(std::vector<float>& audio, float target, float fastAttack, float slowDecay);

// 6) Décimation FIR simple
std::vector<float> simpleFIRDecimate(const std::vector<float>& in, int decim, float cutoffRel = 0.45f);

// 7) Filtre audio final (Biquad Bandpass)
void biquadInitBandpass(Biquad& f, float fs, float f0, float Q);
void biquadProcess(Biquad& f, std::vector<float>& x);
void biquadInitHighpass(Biquad& f, float fs, float f0, float Q);

void transientBoost(std::vector<float>& x, float coeff);

// DETECT PULS - TO DO

// 8) Conversion Float -> PCM16
std::vector<int16_t> floatToPCM(const std::vector<float>& in, float gain);

// 11) demod pipeline
void processSSB_opt(std::vector<std::complex<float>> iq, uint32_t sampleRate,
                    bool upperSideband, std::vector<int16_t> &pcmOut, bool &pulse, int mode);

// 11b) Wide FM demod pipeline — separate from processSSB_opt (rather than
// branching inside it) so the existing, well-tested SSB/pulse-detection path
// used by RDF fox-hunt tracking is untouched. For two-way-radio-style "Wide
// FM" (~25-35kHz channel, ~5kHz deviation), not broadcast FM (200kHz, stereo,
// pre-emphasis) — see soundMode for the volume-mode gain scaling, same
// convention as processSSB_opt.
void processFM_opt(std::vector<std::complex<float>> iq, uint32_t sampleRate,
                    std::vector<int16_t> &pcmOut, int mode);

// Optional: detect if a pulse exists in the last processed block
bool lastBlockContainsPulse();

#endif //SSB_DEMOD_OPT_H
