//
// Created by florent on 27/10/2025.
//

#include "../include/ssb_demod.h"
#include <cmath>
#include <fstream>
#include <algorithm>
#include <cstring>

// ------------------------------------------------------------
// 1. Conversion IQ
// ------------------------------------------------------------
void convertIQ(const unsigned char* buffer, int len, std::vector<std::complex<float>>& iq) {
    int sampCount = len / 2;
    iq.resize(sampCount);
    for (int i = 0; i < sampCount; ++i) {
        float I = (buffer[2 * i]     - 127.4f) / 128.0f;
        float Q = (buffer[2 * i + 1] - 127.4f) / 128.0f;
        iq[i] = std::complex<float>(I, Q);
    }
}

// ------------------------------------------------------------
// 2. Filtre passe-bas
// ------------------------------------------------------------
std::vector<float> makeLowpass(float cutoff, float fs, int N) {
    std::vector<float> h(N);
    int M = N - 1;
    float fc = cutoff / fs;

    for (int n = 0; n < N; n++) {
        int k = n - M/2;
        float sinc = (k == 0) ? 2*M_PI*fc : std::sin(2*M_PI*fc*k)/(float)k;
        float w = 0.54f - 0.46f * std::cos(2*M_PI*n/M); // Hamming
        h[n] = (sinc / M_PI) * w;
    }

    float sum = 0.0f;
    for (auto v : h) sum += v;
    for (auto &v : h) v /= sum;
    return h;
}

// ------------------------------------------------------------
// 3. Application du filtre FIR
// ------------------------------------------------------------
std::vector<std::complex<float>> applyFIR(const std::vector<std::complex<float>>& in, const std::vector<float>& h) {
    int N = h.size();
    int L = in.size();
    std::vector<std::complex<float>> out(L, {0,0});

    for (int i = 0; i < L; ++i) {
        std::complex<float> acc(0,0);
        for (int k = 0; k < N; ++k) {
            if (i - k < 0) break;
            acc += in[i - k] * h[k];
        }
        out[i] = acc;
    }
    return out;
}

// ------------------------------------------------------------
// 4. Démodulation SSB (product detector)
// ------------------------------------------------------------
std::vector<float> demodSSB(const std::vector<std::complex<float>>& sig, bool upperSideband) {
    std::vector<float> audio(sig.size());
    for (size_t i = 0; i < sig.size(); ++i) {
        if (upperSideband)
            audio[i] = sig[i].real() + sig[i].imag();  // USB
        else
            audio[i] = sig[i].real() - sig[i].imag();  // LSB
    }
    return audio;
}

// ------------------------------------------------------------
// 5. Décimation rapide
// ------------------------------------------------------------
std::vector<float> decimate(const std::vector<float>& in, int factor) {
    std::vector<float> out;
    out.reserve(in.size() / factor);
    for (size_t i = 0; i < in.size(); i += factor)
        out.push_back(in[i]);
    return out;
}

// ------------------------------------------------------------
// 6. Conversion float -> PCM16
// ------------------------------------------------------------
std::vector<int16_t> floatToPCM(const std::vector<float>& in, float gain) {
    std::vector<int16_t> pcm(in.size());
    for (size_t i = 0; i < in.size(); ++i) {
        float v = std::clamp(in[i] * gain, -1.0f, 1.0f);
        pcm[i] = static_cast<int16_t>(v * 32767);
    }
    return pcm;
}

// ------------------------------------------------------------
// 7. Pipeline complet IQ -> PCM
// ------------------------------------------------------------
void processSSB(const unsigned char* buffer, int len, uint32_t sampleRate, bool upperSideband, std::vector<int16_t>& pcmOut, float gain) {
    std::vector<std::complex<float>> iq;
    convertIQ(buffer, len, iq);

    std::vector<float> h = makeLowpass(3000.0f, sampleRate, 257);
    auto filt = applyFIR(iq, h);
    auto audio = demodSSB(filt, upperSideband);

    int decim = std::max(1, static_cast<int>(sampleRate / 48000.0f));
    auto audio48k = decimate(audio, decim);
    pcmOut = floatToPCM(audio48k, gain);
}

// ------------------------------------------------------------
// 8. Écriture WAV 16-bit mono
// ------------------------------------------------------------
bool writeWav(const std::string& filename, const std::vector<int16_t>& pcm, int sampleRate) {
    std::ofstream f(filename, std::ios::binary);
    if (!f) return false;

    int dataSize = pcm.size() * sizeof(int16_t);
    int fileSize = 36 + dataSize;

    f.write("RIFF", 4);
    int32_t chunkSize = fileSize;
    f.write(reinterpret_cast<const char*>(&chunkSize), 4);
    f.write("WAVE", 4);

    // fmt chunk
    f.write("fmt ", 4);
    int32_t subchunk1Size = 16;
    int16_t audioFormat = 1;
    int16_t numChannels = 1;
    int32_t byteRate = sampleRate * numChannels * 2;
    int16_t blockAlign = numChannels * 2;
    int16_t bitsPerSample = 16;
    f.write(reinterpret_cast<const char*>(&subchunk1Size), 4);
    f.write(reinterpret_cast<const char*>(&audioFormat), 2);
    f.write(reinterpret_cast<const char*>(&numChannels), 2);
    f.write(reinterpret_cast<const char*>(&sampleRate), 4);
    f.write(reinterpret_cast<const char*>(&byteRate), 4);
    f.write(reinterpret_cast<const char*>(&blockAlign), 2);
    f.write(reinterpret_cast<const char*>(&bitsPerSample), 2);

    // data chunk
    f.write("data", 4);
    f.write(reinterpret_cast<const char*>(&dataSize), 4);
    f.write(reinterpret_cast<const char*>(pcm.data()), dataSize);

    return true;
}

