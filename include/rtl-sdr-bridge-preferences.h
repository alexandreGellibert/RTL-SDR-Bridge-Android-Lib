//
// Created by alexandre on 4/20/25.
//

#ifndef SDR_BRIDGE_PREFERENCES_H
#define SDR_BRIDGE_PREFERENCES_H


#include <string>
#include <chrono>

class Preferences {
public:
    // Singleton access
    static Preferences& getInstance() {
        static Preferences instance;
        return instance;
    }

    // Initialize with values from Kotlin
    void initialize(
            uint32_t centerFrequency,
            uint32_t sampleRate,
            int samplesPerReading,
            int freqFocusRangeKhz,
            int gain,
            long long refreshFFTMs,
            long long refreshPeakMs,
            long long refreshSignalStrengthMs,
            float ssbGain,
            bool dynamicThreshold,
            bool narrowWindow,
            float wwThresholdLVL2,
            float wwThresholdLVL3
    ) {
        centerFrequency_ = centerFrequency;
        sampleRate_ = sampleRate;
        samplesPerReading_ = samplesPerReading;
        freqFocusRangeKhz_ = freqFocusRangeKhz;
        gain_ = gain;
        refreshFFTMs_ = std::chrono::milliseconds(refreshFFTMs);
        refreshPeakMs_ = std::chrono::milliseconds(refreshPeakMs);
        refreshSignalStrengthMs_ = std::chrono::milliseconds(refreshSignalStrengthMs);
        ssbGain_ = ssbGain;
        dynamicThreshold_ = dynamicThreshold;
        narrowWindow_ = narrowWindow;
        wwThresholdLVL2_ = wwThresholdLVL2;
        wwThresholdLVL3_ = wwThresholdLVL3;
        isPrefsInitialized_ = true;
    }

    // Getters
    uint32_t getCenterFrequency() const { return centerFrequency_; }
    uint32_t getSampleRate() const { return sampleRate_; }
    int getFreqFocusRangeKhz() const { return freqFocusRangeKhz_; }
    int getSamplesPerReading() const { return samplesPerReading_; }
    int getGain() const { return gain_; }
    std::chrono::milliseconds getRefreshGraphMs() const { return refreshFFTMs_; }
    std::chrono::milliseconds getRefreshStrengthMs() const { return refreshPeakMs_; }
    std::chrono::milliseconds getBipMaxLengthMs() const { return refreshSignalStrengthMs_; }
    float getSsbGain() const { return ssbGain_; }
    bool getDynamicThreshold() const { return dynamicThreshold_; }
    bool getNarrowWindow() const { return narrowWindow_; }
    float getWwThresholdLVL2() const { return wwThresholdLVL2_; }
    float getWwThresholdLVL3() const { return wwThresholdLVL3_; }
    bool isInitialized() const { return isPrefsInitialized_; }

    // Setters (for updates after initialization if needed)
    void setCenterFrequency(long long value) { centerFrequency_ = value; }
    void setSampleRate(long long value) { sampleRate_ = value; }
    void setFreqFocusRangeKhz(int value) { freqFocusRangeKhz_ = value; }
    void setSamplesPerReading(int value) { samplesPerReading_ = value; }
    void setGain(int value) { gain_ = value; }
    void setRefreshFFTMs(long long value) { refreshFFTMs_ = std::chrono::milliseconds(value); }
    void setRefreshPeakMs(long long value) { refreshPeakMs_ = std::chrono::milliseconds(value); }
    void setRefreshSignalStrengthMs(long long value) { refreshSignalStrengthMs_ = std::chrono::milliseconds(value); }
    void setSsbGain(float value) { ssbGain_ = value; }
    void setDynamicThreshold(bool value) { dynamicThreshold_ = value; }
    void setNarrowWindow(bool value) { narrowWindow_ = value; }
    void setWwThresholdLVL2(float value) { wwThresholdLVL2_ = value; }
    void setWwThresholdLVL3(float value) { wwThresholdLVL3_ = value; }

private:
    Preferences() = default; // Private constructor for singleton

    // Member variables
    uint32_t centerFrequency_ = 0;
    uint32_t sampleRate_ = 0;
    int freqFocusRangeKhz_ = 0;
    int samplesPerReading_ = 0;
    int gain_ = 0;
    std::chrono::milliseconds refreshFFTMs_ = std::chrono::milliseconds(0);
    std::chrono::milliseconds refreshPeakMs_ = std::chrono::milliseconds(0);
    std::chrono::milliseconds refreshSignalStrengthMs_ = std::chrono::milliseconds(0);
    float ssbGain_ = 1.2f;
    bool dynamicThreshold_ = false;
    bool narrowWindow_ = false;
    float wwThresholdLVL2_ = 0.0f;
    float wwThresholdLVL3_ = 0.0f;
    bool isPrefsInitialized_ = false;
};


#endif //SDR_BRIDGE_PREFERENCES_H
