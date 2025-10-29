//
// Created by florent on 27/10/2025.
//

#ifndef SSB_DEMOD_H
#define SSB_DEMOD_H

#include <vector>
#include <complex>
#include <cstdint>

/// Convertit un buffer I/Q brut (uint8 interleaved) en échantillons complexes float
void convertIQ(const unsigned char* buffer, int len, std::vector<std::complex<float>>& iq);

/// Construit un filtre passe-bas FIR simple
std::vector<float> makeLowpass(float cutoff, float fs, int N = 257);

/// Applique un filtre FIR (coefficients réels) sur un signal complexe
std::vector<std::complex<float>> applyFIR(const std::vector<std::complex<float>>& in, const std::vector<float>& h);

/// Démodule un signal SSB (USB ou LSB)
std::vector<float> demodSSB(const std::vector<std::complex<float>>& sig, bool upperSideband = true);

/// Décime le signal par un facteur entier (sans anti-alias, à éviter en production)
std::vector<float> decimate(const std::vector<float>& in, int factor);

/// Convertit un signal float (-1..1) en PCM 16 bits
std::vector<int16_t> floatToPCM(const std::vector<float>& in, float gain = 0.8f);

/// Chaîne complète de traitement SSB IQ -> PCM
void processSSB(const unsigned char* buffer, int len, uint32_t sampleRate, bool upperSideband, std::vector<int16_t>& pcmOut, float gain = 5.0f);

/// Écrit un buffer PCM en fichier WAV 16-bit mono
bool writeWav(const std::string& filename, const std::vector<int16_t>& pcm, int sampleRate = 48000);

#endif // SSB_DEMOD_H
