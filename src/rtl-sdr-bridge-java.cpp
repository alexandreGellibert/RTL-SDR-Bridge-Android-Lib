#include <cmath>
#include <jni.h>
#include <string.h>
#include <__exception/exception.h>
#include <android/log.h>
#include <rtl-sdr-android.h>
#include <fftw3.h>
#include <cstdlib>
#include <__algorithm/max.h>
#include <vector>
#include <chrono>
#include <iostream>
#include <numeric>
#include <complex>
#include <rtl-sdr.h>
#include "../include/rtl-sdr-bridge-preferences.h"
#include <algorithm>
#include <stdexcept>
#include <cstring>

#define LOG_TAG "RTL-SDR Bridge"
#define BASENAME(file) (strrchr(file, '/') ? strrchr(file, '/') + 1 : (strrchr(file, '\\') ? strrchr(file, '\\') + 1 : file))

// Cache the JavaVM globally
JavaVM *gJavaVM = nullptr;

// Get the JavaVM (implementation depends on your setup)
JavaVM *getJavaVM() {
    return gJavaVM;
}

// Function to send logs to Java/Kotlin
void logToJava(const char *message) {
    if (gJavaVM == nullptr) {
        __android_log_print(ANDROID_LOG_ERROR, LOG_TAG, "JavaVM is NULL. Cannot logToJava");
        return;
    }

    JNIEnv *env = nullptr;
    bool attached = false;

    // Check if the thread is already attached to the JVM
    jint res = gJavaVM->GetEnv((void **) &env, JNI_VERSION_1_6);
    if (res == JNI_EDETACHED) {
        // Thread is not attached, so attach it
        if (gJavaVM->AttachCurrentThread(&env, nullptr) == JNI_OK) {
            attached = true; // Mark that we attached the thread
        } else {
            __android_log_print(ANDROID_LOG_ERROR, LOG_TAG, "Failed to attach thread to JVM");
            return;
        }
    } else if (res != JNI_OK) {
        __android_log_print(ANDROID_LOG_ERROR, LOG_TAG, "Failed to get JNIEnv");
        return;
    }

    jclass wrapperClass = env->FindClass("fr/intuite/rtlsdrbridge/RtlSdrBridgeWrapper");
    jmethodID logMethod = env->GetStaticMethodID(wrapperClass, "logFromNative",
                                                 "(Ljava/lang/String;)V");
    jstring jMessage = env->NewStringUTF(message);

    if (jMessage != nullptr) {
        env->CallStaticVoidMethod(wrapperClass, logMethod, jMessage);
        if (env->ExceptionOccurred()) {
            __android_log_print(ANDROID_LOG_ERROR, LOG_TAG,
                                "Exception while calling Java Logs from C++");
            env->ExceptionDescribe();
            env->ExceptionClear();
        }
        env->DeleteLocalRef(jMessage);
    } else {
        __android_log_print(ANDROID_LOG_ERROR, LOG_TAG, "Failed to create Java string");
    }
    if (attached) {
        gJavaVM->DetachCurrentThread();
    }

}

#define LOGD(fmt, ...) do { \
    char buf[1024]; \
    snprintf(buf, sizeof(buf), "[%s:%d] " fmt, BASENAME(__FILE__), __LINE__, ##__VA_ARGS__); \
    if (strlen(buf) >= sizeof(buf) - 1) { \
        __android_log_print(ANDROID_LOG_WARN, LOG_TAG, "Log truncated"); \
    } \
    logToJava(buf); \
} while(0)

// Cache JavaVM during JNI initialization
extern "C" JNIEXPORT jint JNICALL JNI_OnLoad(JavaVM *vm, void *reserved) {
    gJavaVM = vm; // Store the JavaVM pointer
    __android_log_print(ANDROID_LOG_DEBUG, LOG_TAG, "JNI_OnLoad called, JavaVM cached");
    return JNI_VERSION_1_6;
}

static rtlsdr_dev_t *dev = nullptr;
static jobject fftCallbackObj = nullptr; // Global reference to the callback
static jmethodID fftCallbackMethod = nullptr; // Method ID for the callback
static jobject strengthCallbackObj = nullptr; // Global reference to the callback
static jmethodID strengthCallbackMethod = nullptr; // Method ID for the callback
static jobject peakCallbackObj = nullptr; // Global reference to the callback
static jmethodID peakCallbackMethod = nullptr; // Method ID for the callback
static jobject peakNormalizedCallbackObj = nullptr; // Global reference to the callback
static jmethodID peakNormalizedCallbackMethod = nullptr; // Method ID for the callback
static jobject peakFrequencyCallbackObj = nullptr; // Global reference to the callback
static jmethodID peakFrequencyCallbackMethod = nullptr; // Method ID for the callback
static bool isRunning = false;
static bool isCenterFrequencyChanged =false ;

extern "C" JNIEXPORT jboolean JNICALL
Java_fr_intuite_rtlsdrbridge_RtlSdrBridgeWrapper_nativeInitRTL(JNIEnv *env, jobject obj, jint fd, jstring usbfsPath) {

    Preferences &prefs = Preferences::getInstance();
    if (!prefs.isInitialized()) {
        LOGD("Preferences and Settings not yet initialized");
        return false;
    }

    const char *usbfsPathCStr = env->GetStringUTFChars(usbfsPath, nullptr);
    int ret = rtlsdr_open2(&dev, fd, usbfsPathCStr);
    if (ret < 0) {
        LOGD("rtlsdr_open2 failed: %d", ret);
        env->ReleaseStringUTFChars(usbfsPath, usbfsPathCStr);
        return false;
    }

    rtlsdr_dev_t *device = dev;
    int r = 0;
    if (prefs.getCenterFrequency() == 0 ||
        (r = rtlsdr_set_center_freq(device, prefs.getCenterFrequency()) < 0)) {
        LOGD("ERROR: Failed to frequency to %ld with result %d", prefs.getCenterFrequency(), r);
        return false;
    } else {
        LOGD("Init frequency to %ld", prefs.getCenterFrequency());
    }

    if (prefs.getSampleRate() == 0 ||
        (r = rtlsdr_set_sample_rate(device, prefs.getSampleRate())) < 0) {
        LOGD("ERROR: Failed to set sample rate to %lld with result %d", prefs.getSampleRate(), r);
        return false;
    } else {
        LOGD("Init sampling rate to %lld", prefs.getSampleRate());
    }

    if (prefs.getGain() < 0) {
        if (r = rtlsdr_set_tuner_gain_mode(device, 0) < 0) {
            LOGD("WARNING: Failed to enable automatic gain with result %d", r);
        }
    } else {
        /* Enable manual gain */
        if (r = rtlsdr_set_tuner_gain_mode(device, 1) < 0)
            LOGD("WARNING: Failed to enable manual gain with result %d", r);
        if (r = rtlsdr_set_tuner_gain(device, prefs.getGain()) < 0)
            LOGD("WARNING: Failed to set tuner gain with result %d", r);
        else
            LOGD("Tuner gain set to %f dB", prefs.getGain() / 10.0);
    }

//    Unused Parameters
//    rtlsdr_set_freq_correction(_this->openDev, _this->ppm);
//    rtlsdr_set_tuner_bandwidth(_this->openDev, 0);
//    rtlsdr_set_direct_sampling(_this->openDev, _this->directSamplingMode);
//    rtlsdr_set_bias_tee(_this->openDev, _this->biasT);
//    rtlsdr_set_agc_mode(_this->openDev, _this->rtlAgc);
//    rtlsdr_set_tuner_gain(_this->openDev, _this->gainList[_this->gainId]);

    if (r = rtlsdr_reset_buffer(device) < 0)
        LOGD("WARNING: Failed to reset buffers with result %d", r);

    env->ReleaseStringUTFChars(usbfsPath, usbfsPathCStr);
    return true;
}

// Global or static jfloatArray (ensure thread safety if needed)
static jfloatArray result = nullptr;

// Declaration of signal treatment variables
using namespace std::chrono;
// FFT and visiualization
static int integrationCount = 0;
static int integrationPeriod = 8; // Number of packets to integrate in temporal
static int bufferIndex = 0;
static fftwf_complex *circularBuffer = nullptr ;
static int currentSampCount = 0;
//Frequencey of interest management
static float frequence_of_interest = 432000000.0f;
static int index_f_interest = 5;
static float trackingFrequency = 0.0f ;
static float lowTrackingFrequency = - 3000.0f ; //inf tol -3khz
static float upperTrackingFrequency = 1000.0f ; //sup tol +1khz
static std::vector<float> maxPeakAndFrequency ;
static std::chrono::steady_clock::time_point timeOfLastMaxPeak;
static std::chrono::steady_clock::time_point timeOfLastMaxPeakUpdate;
//Peak variables
static float refMagnitude = 50000.0f;
static float refPower = refMagnitude * refMagnitude;
//Signal strength variables
std::vector<time_point<high_resolution_clock>> signalTimeTable;
static std::vector<int> signalStrengthIndexBuffer ;
static int signalStrengthIndexRemanance = 30 ;
static int indexsignalStrengthIndexBuffer = 0 ;
static int signalWeak = 1 ;
static int signalMedium = 5 ;
static int signalStrong = 20 ;
static int remananceWeak = 0 ;
static int remananceMedium = 0 ;
static int remananceStrong = 0 ;
//Noise evaluation
static std::vector<float> noisePercentileBuffer ;
static std::vector<float> noiseMedianBuffer ;
static int noiseBufferSize = 15 ;
static int indexNoiseBuffer = 0 ;
//Level 1 variables
static std::vector<float> peakBuffer ;
static std::vector<float> peakNormalizedBuffer ;
static int peakRemanance = 50 ;
static int indexPeakBuffer = 0 ;
//Level 2 variables
static std::vector<float> circularPowerDBnormalized ; // Used to stock power DB (not integrated) centered on 0
static int indexCircularPowerDB = 0; // Circulate within circularPower
const int integrationPowerPeriod = 21; // Number of loop tocked for power DB normalized work on !
static int longueurTrace = 6 ; // Must be < at integrationPower period : circular
static int previousIndexLvl2 = 666 ;
static int lvl2_repet = 0 ;
static int lvl2_repet_threshold =3 ;
static float maxlvl2 = 0.0f ;
//Level 3 variables
int longueurTrace_zigzag = 13 ; //Must be < at integrationPower period : circular
static int previousIndexLvl3 = 666 ;
static int lvl3_repet = 0 ;
static int lvl3_repet_threshold =5 ;
static float maxlvl3 = 0.0f ;


extern "C" JNIEXPORT void JNICALL
Java_fr_intuite_rtlsdrbridge_RtlSdrBridgeWrapper_nativeReadAsync(
        JNIEnv *env, jobject thiz,
        jobject fftCallback,
        jobject signalStrengthCallback,
        jobject peakCallback,
        jobject peakNormalizedCallback,
        jobject peakFrequencyCallback) {

    if (dev == nullptr) {
        LOGD("Device not initialized");
        return;
    }

    // Store the fftCallback globally
    fftCallbackObj = env->NewGlobalRef(fftCallback);
    jclass fftCallbackClass = env->GetObjectClass(fftCallback);
    fftCallbackMethod = env->GetMethodID(fftCallbackClass, "invoke", "([F)V");

    strengthCallbackObj = env->NewGlobalRef(signalStrengthCallback);
    jclass strengthCallbackClass = env->GetObjectClass(signalStrengthCallback);
    strengthCallbackMethod = env->GetMethodID(strengthCallbackClass, "invoke", "(I)V");

    peakCallbackObj = env->NewGlobalRef(peakCallback);
    jclass peakCallbackClass = env->GetObjectClass(peakCallback);
    peakCallbackMethod = env->GetMethodID(peakCallbackClass, "invoke","(F)V");

    peakNormalizedCallbackObj = env->NewGlobalRef(peakNormalizedCallback);
    jclass peakNormalizedCallbackClass = env->GetObjectClass(peakNormalizedCallback);
    peakNormalizedCallbackMethod = env->GetMethodID(peakNormalizedCallbackClass, "invoke","(F)V");

    peakFrequencyCallbackObj = env->NewGlobalRef(peakFrequencyCallback);
    jclass peakFrequencyCallbackClass = env->GetObjectClass(peakFrequencyCallback);
    peakFrequencyCallbackMethod = env->GetMethodID(peakFrequencyCallbackClass, "invoke","(J)V");

    // Callback for rtlsdr_read_async
    auto rtlsdrCallback = [](unsigned char *buffer, uint32_t len, void *ctx) {

        if (!isRunning) {
            return;
        }

        JNIEnv *env;
        // Check if the thread is already attached
        jint result2 = getJavaVM()->GetEnv((void **) &env, JNI_VERSION_1_6);

        uint32_t centerFrequency = Preferences::getInstance().getCenterFrequency();
        uint32_t sampleRate = Preferences::getInstance().getSampleRate();

        int sampCount = len / 2;

        fftwf_complex *signal = fftwf_alloc_complex(sampCount);
        fftwf_complex *fft_signal = fftwf_alloc_complex(sampCount);

        if (currentSampCount != sampCount) {
            fftwf_free(circularBuffer);
            currentSampCount = sampCount;
            circularBuffer = fftwf_alloc_complex(currentSampCount * integrationPeriod);
        }


        for (int i = 0; i < sampCount && i * 2 < len; i++) {
            signal[i][0] = (float) (buffer[2 * i] - 127.4) / 128.0;
            signal[i][1] = (float) (buffer[2 * i + 1] - 127.4) / 128.0;

            // Update circular buffer
            circularBuffer[(bufferIndex * sampCount) + i][0] = signal[i][0];
            circularBuffer[(bufferIndex * sampCount) + i][1] = signal[i][1];
        }

        bufferIndex = (bufferIndex + 1) % integrationPeriod;
        integrationCount++;

        // Integrate the buffer
        fftwf_complex *signal_integrated = fftwf_alloc_complex(sampCount);
        fftwf_complex *fft_signal_integrated = fftwf_alloc_complex(sampCount);

        // Remettre le tableau à zéro
        for (int i = 0; i < sampCount; ++i) {
            signal_integrated[i][0] = 0.0f; // Partie réelle
            signal_integrated[i][1] = 0.0f; // Partie imaginaire
        }

        // Perform integration
        //TO DO : integrate on timelapse and not on number of loop : if perf is different chat happens ?
        if (integrationCount >= integrationPeriod) {
            for (int i = 0; i < sampCount; i++) {
                for (int j = 0; j < integrationPeriod; j++) {
                    signal_integrated[i][0] += circularBuffer[j * sampCount + i][0];
                    signal_integrated[i][1] += circularBuffer[j * sampCount + i][1];
                }
                //LOGD("passe dans la boucle integration");
                signal_integrated[i][0] /= integrationPeriod;
                signal_integrated[i][1] /= integrationPeriod;
            }
        }

        // Perform FFT on signal
        fftwf_plan plan1 = fftwf_plan_dft_1d(sampCount, signal, fft_signal, FFTW_FORWARD,
                                             FFTW_ESTIMATE);
        fftwf_execute(plan1);
        // AND clean up
        fftwf_destroy_plan(plan1);


        // Perform FFT on signal integrated
        fftwf_plan plan2 = fftwf_plan_dft_1d(sampCount, signal_integrated, fft_signal_integrated,
                                             FFTW_FORWARD, FFTW_ESTIMATE);
        fftwf_execute(plan2);
        // AND clean up
        fftwf_destroy_plan(plan2);

        //Define power inst and power integrated
        float power[sampCount];
        float power_shifted[sampCount];
        float power_integrated[sampCount];
        float power_integrated_shifted[sampCount];

        //Compute power spectrum directly on signal integrated: power = real^2 + imag^2
        for (int i = 0; i < sampCount; i++) {
            power[i] = fft_signal[i][0] * fft_signal[i][0] + fft_signal[i][1] * fft_signal[i][1];
        }
        for (int i = 0; i < sampCount; i++) {
            power_integrated[i] = fft_signal_integrated[i][0] * fft_signal_integrated[i][0] +
                                  fft_signal_integrated[i][1] * fft_signal_integrated[i][1];
        }
        //FREE
        fftwf_free(signal);
        fftwf_free(signal_integrated);
        fftwf_free(fft_signal);
        fftwf_free(fft_signal_integrated);

        // Shift the power spectrum to center it (FFT shift)
        int half = sampCount / 2;
        for (int i = 0; i < half; i++) {
            power_shifted[i] = power[i + half]; // Negative frequencies
            power_shifted[i + half] = power[i]; // Positive frequencies
        }
        for (int i = 0; i < half; i++) {
            power_integrated_shifted[i] = power_integrated[i + half]; // Negative frequencies
            power_integrated_shifted[i + half] = power_integrated[i]; // Positive frequencies
        }

        if (result == nullptr) {
            result = env->NewFloatArray(sampCount);
            // Make it a global reference to persist across calls
            result = (jfloatArray) env->NewGlobalRef(result);
        }

        //NEXT : give power in db
        env->SetFloatArrayRegion(result, 0, sampCount, power_integrated_shifted);
        // Call the Kotlin fftCallback
        env->CallVoidMethod(fftCallbackObj, fftCallbackMethod, result);

        // BANDWITH for PEAK EVALUATION
        // Compute signal strength in subrange = wide window = WW(fc ± 10 kHz)
        float freqPerBin = static_cast<float>(sampleRate) / sampCount;  // Hz per bin
        float lowerWWBound = centerFrequency - Preferences::getInstance().getFreqFocusRangeKhz() *
                                             1000.0f;  // fc - 10 kHz
        float upperWWBound = centerFrequency + Preferences::getInstance().getFreqFocusRangeKhz() *
                                             1000.0f;  // fc + 10 kHz

        // Find bin indices for the subrange (after shifting)
        int lowerWWIndex = static_cast<int>((lowerWWBound - (centerFrequency - sampleRate / 2)) /
                                          freqPerBin);
        int upperWWIndex = static_cast<int>((upperWWBound - (centerFrequency - sampleRate / 2)) /
                                          freqPerBin);
        lowerWWIndex = std::max(0, lowerWWIndex);
        upperWWIndex = sampCount - 1 > upperWWIndex ? upperWWIndex : sampCount - 1;
        int WW_size = upperWWIndex - lowerWWIndex + 1;

        // SIGNAL PEAK MEASURE
        // Look for signal peak (in case of signal strong enough)
        float totalPower = 0.0f;
        float total_integrated_Power = 0.0f;
        float dB = -130.0f; //temp value for max research
        float peakDb = -130.0f; //peak on integrated signal
        float peakNormalized = 0.0f ; //peak normalized
        int index_peak = 0;
        float WW_average_power = -130.0f; //TO DO remonter en static si on veut moyenner un peu ???
        float WW_average_integrated_power = -130.0f;

        float power_integrated_shifted_DB[WW_size];

        for (int i = lowerWWIndex; i <= upperWWIndex; i++) {
            dB = 10 * log10(power_integrated_shifted[i] / refPower);
            if (dB > peakDb) {
                peakDb = dB;
                index_peak = i - lowerWWIndex;
            }

            totalPower += 10 * log10(power_shifted[i] / refPower);
            power_integrated_shifted_DB[i - lowerWWIndex] =
                    10 * log10(power_integrated_shifted[i] / refPower);
            total_integrated_Power += 10 * log10(power_integrated_shifted[i] / refPower);
        }
        WW_average_power = (totalPower / WW_size); // average defined on signal not integrated
        WW_average_integrated_power = (total_integrated_Power / WW_size);

        //Prepare peak buffer for remanance.
        if (peakBuffer.empty()) {
            peakBuffer.assign(peakRemanance, -130.0);
        }
        //IDEM prepare peak normalized for remanance
        if (peakNormalizedBuffer.empty()) {
            peakNormalizedBuffer.assign(peakRemanance, 0.0);
        }
        //Stock peakDb.
        peakBuffer[indexPeakBuffer] = peakDb;

        //LEVEL 1 SET-UP: peak reading as something appearing on top of noise
        //Noise statistics calculation
        //Init
        if(noisePercentileBuffer.empty()){
            noisePercentileBuffer.assign(noiseBufferSize,0.0);
        }
        if(noiseMedianBuffer.empty()){
            noiseMedianBuffer.assign(noiseBufferSize,0.0);
        }

        //Calculation !
        //Do not take all the data to have some out of peak values
        if (integrationCount % 5 == 0) {
            float pourcentage = 2.0f; // Par exemple, 10%
            int indexOfPercentile = static_cast<int>(std::floor(static_cast<double>(WW_size) * (1 - pourcentage / 100)));
            int indexOfMedian = static_cast<int>(std::floor(static_cast<double>(WW_size) * 0.5)); //Middle ;)
            //LOGD("index = %d", index) ;
            std::sort(power_integrated_shifted_DB, power_integrated_shifted_DB + WW_size);
            noisePercentileBuffer[indexNoiseBuffer] = power_integrated_shifted_DB[indexOfPercentile] ;
            noiseMedianBuffer[indexNoiseBuffer] = power_integrated_shifted_DB[indexOfMedian] ;
            //LOGD("valeur percentilebuffer = %f", percentileBuffer[indexPercentileBuffer]);
            indexNoiseBuffer = (indexNoiseBuffer + 1) % noiseBufferSize;
        }
        //Define temp tables to sort data and get percentile and median
        float percentileBufferSorted [noiseBufferSize] ;
        float medianBufferSorted [noiseBufferSize] ;
        for (int k = 0 ; k<noiseBufferSize;k++){
            percentileBufferSorted[k]=noisePercentileBuffer[k];
        }
        for (int k = 0 ; k<noiseBufferSize;k++){
            medianBufferSorted[k]=noiseMedianBuffer[k];
        }
        std::sort(percentileBufferSorted, percentileBufferSorted + noiseBufferSize);
        std::sort(medianBufferSorted, medianBufferSorted + noiseBufferSize);
        //Calculate data
        float noisePercentile = percentileBufferSorted[noiseBufferSize-static_cast<int>(std::floor(static_cast<double>(noiseBufferSize)/2))];
        float noiseMedian = medianBufferSorted[noiseBufferSize-static_cast<int>(std::floor(static_cast<double>(noiseBufferSize)/2))];
        float noiseSigma = noisePercentile - noiseMedian ; //not really sigma but...

        peakNormalized = peakDb - noiseMedian ;
        //IDEM Stock peakNormalized
        peakNormalizedBuffer[indexPeakBuffer] = peakNormalized ;
        //increment index
        indexPeakBuffer = (indexPeakBuffer + 1) % peakRemanance;

        //CHECK
        if (integrationCount % 500 == 0) {
            LOGD("Every 500 round. gap = %f", noiseSigma);
        }

        //INIT of tables for narrow window

        if(trackingFrequency == 0.0f){
            trackingFrequency = centerFrequency ;
        }
        if(isCenterFrequencyChanged){
            trackingFrequency = centerFrequency ;
            isCenterFrequencyChanged = false ;
        }
        if (maxPeakAndFrequency.empty()){
            maxPeakAndFrequency = {-130.0f, static_cast<float>(centerFrequency)};
        }

        //LEVEL 2 SET-UP : set up. trace des instantannés : ressemble à multiplier les signaux ;). Bon paramètre : integration period = 13, moyenne 5.5. Sur instantanné. Sur longueur 6 avec déclencheur à 7 = ok aussi ! Mettre de la récurrence.
        //Define circular power size
        circularPowerDBnormalized.resize(WW_size * integrationPowerPeriod);

        //FEED circular power matrix with lines centered on 0 as average of power en instantaneous signal
        for (int i = 0; i < WW_size && i * 2 < len; i++) {
            int indexGlobal = i + lowerWWIndex;
            // Update circular buffer with power
            circularPowerDBnormalized[(indexCircularPowerDB * WW_size) + i] = 10 * log10(power_shifted[indexGlobal] / refPower) - WW_average_power;
        }
        //Calculate trace
        std::vector<float> trace(WW_size, 0.0f);
        float maxTrace = 0.0f;
        int index_lvl2 = 0;
        for (int i = 0; i < WW_size; i++) {
            for (int k = 0; k < longueurTrace; k++) {
                trace[i] += circularPowerDBnormalized[i + ((indexCircularPowerDB - k+integrationPowerPeriod)%integrationPowerPeriod) * WW_size];
            }
            if (trace[i] > maxTrace) {
                maxTrace = trace[i];
                index_lvl2 = i;
            }
        }
        //Used only to try to define proper value
        if (maxTrace>maxlvl2){
            maxlvl2 = maxTrace ;
        }

        // LEVEL 3 SET-UP : éclairs ou cones. toper les obliques. A finir de régler !!!
        std::vector<float> trace_zigzag(WW_size, 0.0f);
        float maxTrace_zigzag = 0.0f;
        int index_lvl3 = 0;
        for (int i = 0; i < WW_size; i++) {
            int i_zigzag = i ;
            float max_zigzag = -10.0f ;
            for (int k = 0; k < longueurTrace_zigzag; k++) {
                if (k==0){
                    max_zigzag=circularPowerDBnormalized[i + ((indexCircularPowerDB - k + integrationPowerPeriod) %integrationPowerPeriod) *WW_size];
                }
                else{
                    for (int t = std::max(i_zigzag-1, 0);t<=std::min(i_zigzag+1,longueurTrace_zigzag);t++) {
                         if(max_zigzag<circularPowerDBnormalized[t + ((indexCircularPowerDB - k + integrationPowerPeriod) %integrationPowerPeriod) *WW_size]){
                             max_zigzag = circularPowerDBnormalized[t + ((indexCircularPowerDB - k + integrationPowerPeriod) %integrationPowerPeriod) *WW_size];
                             i_zigzag = t ;
                         }

                    }
                }
                trace_zigzag[i] += max_zigzag ;
            }
            if (trace_zigzag[i] > maxTrace_zigzag) {
                maxTrace_zigzag = trace_zigzag[i];
                index_lvl3 = i;
            }
        }
        //Used only to try to define proper value for lvl3
        if (maxTrace_zigzag>maxlvl3){
            maxlvl3 = maxTrace_zigzag ;
        }

        // TO DO Level 4 : stat sur répartitions

        // signalStrengthIndex get value 0 : no signal ; 1 ; weak signal may be no signal, remanance short ; medium signal, remanance normal ; strong signal, remanance normal.
        //CALCULATION of Signal Strength
        int signalEval = 0 ;

        //LEVEL 1 CALCULATION : first
        //First condition is a bit arbitral. To try link with noise "form"
        if (peakDb > noiseMedian + 2.5 * noiseSigma) {
            if (peakDb>maxPeakAndFrequency[0]) {
                maxPeakAndFrequency[0] = peakDb ;
                maxPeakAndFrequency[1] = index_peak * freqPerBin +lowerWWBound; //définir les règles de lautocalibration !!!
                timeOfLastMaxPeak = std::chrono::steady_clock::now();
                LOGD("Level 1 déclenche. frequence = %f", trackingFrequency);
            }
            signalEval = signalStrong + 1 ;
            LOGD("Level 1 déclenche. peak = %f", peakDb);
            LOGD("Level 1 déclenche. valeur noiseMedian = %f", noiseMedian);
            LOGD("Level 1 déclenche. valeur noiseSigma = %f", noiseSigma);
        }

        auto maintenant = std::chrono::steady_clock::now();
        auto delaiDepuisPeak = std::chrono::duration_cast<std::chrono::milliseconds>(maintenant - timeOfLastMaxPeak);
        //auto delai2 = std::chrono::duration_cast<std::chrono::milliseconds>(timeOfLastMaxPeak - timeOfLastMaxPeakUpdate);
        auto delaiDepuisPeak_ms = delaiDepuisPeak.count() ;

        //PB : déclenche alors que pas de peak !!!
        if (timeOfLastMaxPeakUpdate < timeOfLastMaxPeak && delaiDepuisPeak_ms > 300){
            trackingFrequency = maxPeakAndFrequency[1];
            timeOfLastMaxPeakUpdate = std::chrono::steady_clock::now();
            maxPeakAndFrequency[0] = -130.0f ;
            LOGD("MAJ frequence = %f", trackingFrequency);
        }

        //LEVEL 2 CALCULATION :
        //Look for narrow window if narrow window and if nothing look large
        //init
        if(previousIndexLvl2==666){
            previousIndexLvl2 = index_lvl2 ;
        }

        //dynamic threshold
        if(Preferences::getInstance().getDynamicThreshold()){
            if (noiseSigma>1) {
                Preferences::getInstance().setWwThresholdLVL2(1.5 * noiseSigma);
            }
        }

        float wwThresholdLVL2 = Preferences::getInstance().getWwThresholdLVL2();
        //evaluation of max trace over trigger repetition
        if(maxTrace > wwThresholdLVL2 * longueurTrace && std::abs(previousIndexLvl2-index_lvl2)<4){
            lvl2_repet++;
        }
        else{
            if (Preferences::getInstance().getNarrowWindow()){
                float NWthresholdLVL2 = 0.95 * wwThresholdLVL2; // Peut-être un peu offensif mais détecte super bien
                float freqLvl2 = index_lvl2*freqPerBin+lowerWWBound;

                if (maxTrace > NWthresholdLVL2 * longueurTrace && std::abs(previousIndexLvl2-index_lvl2)<4 &&lowTrackingFrequency+trackingFrequency<freqLvl2 && freqLvl2<upperTrackingFrequency + trackingFrequency){
                    lvl2_repet++;
                }
                else {
                    lvl2_repet = 0;
                    previousIndexLvl2 = 666;
                }
            }
            else {
                lvl2_repet = 0;
                previousIndexLvl2 = 666;
            }
        }

        //Signal eval calculation with lvl2 repetition
        if(lvl2_repet>lvl2_repet_threshold){
            //remanancelvl2 = 20 ;
            LOGD("LVL2 déclenche. signal lvl2 = %f", maxTrace/longueurTrace);
            LOGD("LVL2 déclenche. gap = %f", noiseSigma);
            //LOGD("Index lvl2: %d", index_lvl2);
            //LOGD("lvl2_repet: %d", lvl2_repet);
            signalEval = signalEval+lvl2_repet-lvl2_repet_threshold;
        }

        //LEVEL 3 CALCULATION :
        //init
        if(previousIndexLvl3==666){
            previousIndexLvl3 = index_lvl3 ;
        }

        //dynamic threshold
        if(Preferences::getInstance().getDynamicThreshold()){
            if(noiseSigma>1) {
                Preferences::getInstance().setWwThresholdLVL3(1.8 * noiseSigma);
            }
        }

        //11 pas mal mais déclanche un poil trop. autre idée allonger la trace

        float wwThresholdLVL3 = Preferences::getInstance().getWwThresholdLVL3();
        if (maxTrace_zigzag > wwThresholdLVL3 * longueurTrace_zigzag && std::abs(previousIndexLvl3-index_lvl3)<10) {
            lvl3_repet++;
        }
        else{
            if (Preferences::getInstance().getNarrowWindow()) {
                float NWthresholdLVL3 = 0.95 * wwThresholdLVL3;
                float freqLvl3 = index_lvl3 * freqPerBin + lowerWWBound;
                if (maxTrace_zigzag > wwThresholdLVL3 * longueurTrace_zigzag &&
                    std::abs(previousIndexLvl3 - index_lvl3) < 10 &&
                    lowTrackingFrequency + trackingFrequency < freqLvl3 &&
                    freqLvl3 < upperTrackingFrequency + trackingFrequency) {
                    lvl3_repet++;
                } else {
                    lvl3_repet = 0;
                    previousIndexLvl3 = 666;
                }
            }
            else {
                    lvl3_repet = 0;
                    previousIndexLvl3 = 666;
            }
        }

        //Signal eval calculation with lvl3 repetition
        if(lvl3_repet>lvl3_repet_threshold){ //3 pas mal mais déclanche un poil trop. autre idée allonger la trace
                //remanancelvl3 = 20 ;
                LOGD("LVL3 déclenche. signal lvl3 = %f", maxTrace_zigzag/longueurTrace_zigzag);
                //LOGD("Index lvl3: %d", index_lvl3);
                //LOGD("lvl3_repet: %d", lvl3_repet);
                signalEval=signalEval+lvl3_repet-lvl3_repet_threshold;
        }


        //tableau des occurences des signaux avec valeur et timestamp. weak si repet faible. medium si repet faible et régulier ou si repet fort. fort si repet fort et régulier.
        //A optimiser, pas besoin de passer partout
        if (signalTimeTable.size()==0){
            signalTimeTable.push_back(high_resolution_clock::now());
            signalTimeTable.push_back(high_resolution_clock::now());
        }

        //signal strength calculation
        int signalStrengthIndex = 0 ;

        if (remananceStrong>0){
            signalStrengthIndex = 3 ;
            remananceStrong = std::max(0, remananceStrong-1);
            remananceMedium = std::max(0, remananceMedium-1);
            remananceWeak = std::max(0, remananceWeak-1);
        }
        else {
            // Calculer les durées
            auto duration1 = high_resolution_clock::now() - signalTimeTable[0];
            auto duration2 = signalTimeTable[0] - signalTimeTable[1];

            // Convertir les durées en un type numérique (par exemple, en millisecondes)
            auto duration1_ms = duration_cast<milliseconds>(duration1).count();
            auto duration2_ms = duration_cast<milliseconds>(duration2).count();

            // Calculer la valeur absolue de la différence
            long long abs_diff = std::abs(duration1_ms - duration2_ms);

            if (signalEval >= signalStrong ||
                signalMedium <= signalEval && signalEval < signalStrong && abs_diff < 300) {
                signalStrengthIndex = 3 ;
                remananceStrong = 30;
                if (duration1_ms>666) {
                    signalTimeTable.push_back(high_resolution_clock::now());
                }
            }
            else{
                if(remananceMedium >0){
                    signalStrengthIndex = 2 ;
                    remananceMedium = std::max(0, remananceMedium-1);
                    remananceWeak = std::max(0, remananceWeak-1);
                }
                else{

                    // Calculer les durées
                    auto duration1 = high_resolution_clock::now() - signalTimeTable[0];
                    auto duration2 = signalTimeTable[0] - signalTimeTable[1];

                    // Convertir les durées en un type numérique (par exemple, en millisecondes)
                    auto duration1_ms = duration_cast<milliseconds>(duration1).count();
                    auto duration2_ms = duration_cast<milliseconds>(duration2).count();

                    // Calculer la valeur absolue de la différence
                    long long abs_diff = std::abs(duration1_ms - duration2_ms);

                    //TO DO rajouter duration1_ms > 300 pour pas chercher la confirmation
                    if (signalEval >= signalMedium ||
                        signalWeak <= signalEval && signalEval < signalMedium && abs_diff < 300) {
                        signalStrengthIndex = 2 ;
                        remananceMedium = 15;
                        if (duration1_ms>666) {
                            signalTimeTable.push_back(high_resolution_clock::now());
                        }
                    }
                    else{
                        if(remananceWeak>0){
                            signalStrengthIndex=1 ;
                            remananceWeak = std::max(0, remananceWeak-1);

                        }
                        else{
                            if (signalEval>=signalWeak){
                                remananceWeak =5 ;
                                signalStrengthIndex=1 ;
                                if (duration1_ms>666) {
                                    signalTimeTable.push_back(high_resolution_clock::now());
                                }
                            }
                            else{
                                remananceStrong = std::max(0, remananceStrong-1);
                                remananceMedium = std::max(0, remananceMedium-1);
                                remananceWeak = std::max(0, remananceWeak-1);
                            }

                        }
                    }
                }
            }
        }

        //Initialization of signalStrengthIndex table
        if (signalStrengthIndexBuffer.empty()) {
            signalStrengthIndexBuffer.assign(signalStrengthIndexRemanance, 0);
        }
        //Fulfill of table
        signalStrengthIndexBuffer[indexsignalStrengthIndexBuffer] = signalStrengthIndex ;
        indexsignalStrengthIndexBuffer = (indexsignalStrengthIndexBuffer + 1) % signalStrengthIndexRemanance;

        //Calculate strenghtIndex to be sent
        auto signalStrengthIndexMaxIter = std::max_element(signalStrengthIndexBuffer.begin(), signalStrengthIndexBuffer.end());
        int signalStrengthIndexSent = *signalStrengthIndexMaxIter ;
        //Calculate if we have not quickly some stronger signal before sending "Signal probable"
        int sum = std::accumulate(signalStrengthIndexBuffer.begin(), signalStrengthIndexBuffer.end(), 0);
        double signalStrengthIndexMean = static_cast<double>(sum) / signalStrengthIndexBuffer.size();
        if (*signalStrengthIndexMaxIter  == 1 && signalStrengthIndexMean<0.8){
            int signalStrengthIndexSent = 0 ;
        }

        if (integrationCount%1000==0){
            LOGD("time tous les 1000 tours");
        }

        // Send signal typologyto Java
        env->CallVoidMethod(strengthCallbackObj, strengthCallbackMethod, signalStrengthIndexSent);
        // Stick the peak if goes down
        auto peakRemanantMaxIter = std::max_element(peakBuffer.begin(), peakBuffer.end());
        float peakRemanantMax = *peakRemanantMaxIter ;
        // IDEM for peakNormalized
        auto peakNormalizedMaxIter = std::max_element(peakNormalizedBuffer.begin(), peakNormalizedBuffer.end());
        float peakNormalizedMax = *peakNormalizedMaxIter ;
        // Send peak of signal to Java
        env->CallVoidMethod(peakCallbackObj, peakCallbackMethod, peakRemanantMax);
        // Send normalized peak of signal to Java
        env->CallVoidMethod(peakNormalizedCallbackObj, peakNormalizedCallbackMethod, peakNormalizedMax);
        // Send new frequency tracking to Java
        env->CallVoidMethod(peakFrequencyCallbackObj, peakFrequencyCallbackMethod, static_cast<long>(std::round(trackingFrequency)));

        indexCircularPowerDB = (indexCircularPowerDB + 1) % integrationPowerPeriod;
    };

    rtlsdr_reset_buffer(dev);
    isRunning = true;
    // Start asynchronous reading
    rtlsdr_read_async(dev, rtlsdrCallback, nullptr, 0, Preferences::getInstance().getSamplesPerReading() * 2);
}

// Function to cancel asynchronous reading
extern "C" JNIEXPORT void JNICALL
Java_fr_intuite_rtlsdrbridge_RtlSdrBridgeWrapper_nativeCancelAsync(JNIEnv *env, jobject thiz) {
    if (dev != nullptr) {
        int result = rtlsdr_cancel_async(dev);
        if (result == 0) {
            isRunning = false;
        }
        LOGD("rtlsdr_cancel_async cancel result: %d", result);
    }
}

// Close: Cleanup device
extern "C" JNIEXPORT void JNICALL
Java_fr_intuite_rtlsdrbridge_RtlSdrBridgeWrapper_nativeCloseRTL(JNIEnv *env, jobject obj) {
    isRunning = false;
    if (dev != nullptr) {
        rtlsdr_close(dev);
        dev = nullptr;
    }
    if (fftCallbackObj != nullptr) {
        env->DeleteGlobalRef(fftCallbackObj);
        fftCallbackObj = nullptr;
    }
    if (strengthCallbackObj != nullptr) {
        env->DeleteGlobalRef(strengthCallbackObj);
        strengthCallbackObj = nullptr;
    }
    if (peakFrequencyCallbackObj != nullptr) {
        env->DeleteGlobalRef(peakFrequencyCallbackObj);
        peakFrequencyCallbackObj = nullptr;
    }
    if (peakCallbackObj != nullptr) {
        env->DeleteGlobalRef(peakCallbackObj);
        peakCallbackObj = nullptr;
    }
    if (peakNormalizedCallbackObj != nullptr) {
        env->DeleteGlobalRef(peakNormalizedCallbackObj);
        peakNormalizedCallbackObj = nullptr;
    }
    if (result != nullptr) {
        env->DeleteGlobalRef(result);
        result = nullptr;
    }
    LOGD("Device closed");
}

extern "C" JNIEXPORT void JNICALL
Java_fr_intuite_rtlsdrbridge_RtlSdrBridgeWrapper_nativeSetFrequency(JNIEnv *env, jobject obj,
                                                                 jlong frequency) {
    Preferences::getInstance().setCenterFrequency(frequency);
    if (dev == nullptr) {
        return;
    }

    if (rtlsdr_set_center_freq(dev, (uint32_t) frequency) < 0) {
        LOGD("ERROR: Failed to frequency to %ld", frequency);
    } else {
        LOGD("New Frequency %ld", frequency);
    }
    isCenterFrequencyChanged = true ;
}

extern "C" JNIEXPORT void JNICALL
Java_fr_intuite_rtlsdrbridge_RtlSdrBridgeWrapper_nativeSetSampleRate(JNIEnv *env, jobject obj,
                                                                  jlong sampleRate) {
    Preferences::getInstance().setSampleRate(sampleRate);
    if (dev == nullptr) {
        return;
    }
    int r;
    if ((r = rtlsdr_set_sample_rate(dev, (uint32_t) sampleRate)) < 0) {
        LOGD("ERROR: Failed to set sample rate to %ld with error %d", sampleRate, r);
    } else {
        LOGD("New Sample rate %ld", sampleRate);
    }
}

extern "C" JNIEXPORT void JNICALL
Java_fr_intuite_rtlsdrbridge_RtlSdrBridgeWrapper_nativeSetGain(JNIEnv *env, jobject obj, jint gain) {
    Preferences::getInstance().setGain(gain);
    if (dev == nullptr) {
        return;
    }
    int r;
    if ((r = rtlsdr_set_tuner_gain(dev, (int) gain)) < 0) {
        LOGD("ERROR: Failed to set gain to %d with error %d", gain, r);
    } else {
        LOGD("New Gain %d", gain);
    }
}

extern "C" JNIEXPORT void JNICALL
Java_fr_intuite_rtlsdrbridge_RtlSdrBridgeWrapper_nativeSetSamplesPerReading(JNIEnv *env, jobject obj,
                                                                         jint samplesPerReading) {
    Preferences::getInstance().setSamplesPerReading(samplesPerReading);
    LOGD("New Samples Per Reading %d", samplesPerReading);
}

extern "C" JNIEXPORT void JNICALL
Java_fr_intuite_rtlsdrbridge_RtlSdrBridgeWrapper_nativeSetFrequencyFocusRange(JNIEnv *env, jobject obj,
                                                                           jint frequencyFocusRange) {
    Preferences::getInstance().setSamplesPerReading(frequencyFocusRange);
    LOGD("New Frequency Focus Range %ld", frequencyFocusRange);
}

extern "C" JNIEXPORT void JNICALL
Java_fr_intuite_rtlsdrbridge_RtlSdrBridgeWrapper_nativeSetRefreshFFTMs(JNIEnv *env, jobject obj,
                                                                    jlong refreshFFTMs) {
    Preferences::getInstance().setRefreshFFTMs(refreshFFTMs);
    LOGD("New Refresh FFT period in ms %ld", refreshFFTMs);
}

extern "C" JNIEXPORT void JNICALL
Java_fr_intuite_rtlsdrbridge_RtlSdrBridgeWrapper_nativeSetRefreshPeakMs(JNIEnv *env, jobject obj,
                                                                     jlong refreshPeakMs) {
    Preferences::getInstance().setRefreshPeakMs(refreshPeakMs);
    LOGD("New Refresh Peak period in ms %ld", refreshPeakMs);
}

extern "C" JNIEXPORT void JNICALL
Java_fr_intuite_rtlsdrbridge_RtlSdrBridgeWrapper_nativeSetRefreshSignalStrengthMs(JNIEnv *env,
                                                                               jobject obj,
                                                                               jlong refreshSignalStrengthMs) {
    Preferences::getInstance().setRefreshSignalStrengthMs(refreshSignalStrengthMs);
    LOGD("New Refresh Signal Strength %ld", refreshSignalStrengthMs);
}

extern "C" JNIEXPORT void JNICALL
Java_fr_intuite_rtlsdrbridge_RtlSdrBridgeWrapper_nativeSetIsMuted(JNIEnv *env, jobject obj,
                                                               jboolean isMuted) {
    Preferences::getInstance().setIsMuted(isMuted);
}

extern "C" JNIEXPORT void JNICALL
Java_fr_intuite_rtlsdrbridge_RtlSdrBridgeWrapper_nativeSetDynamicThreshold(JNIEnv *env, jobject obj,
                                                                        jboolean dynamicThreshold) {
    Preferences::getInstance().setDynamicThreshold(dynamicThreshold);
}

extern "C" JNIEXPORT void JNICALL
Java_fr_intuite_rtlsdrbridge_RtlSdrBridgeWrapper_nativeSetNarrowWindow(JNIEnv *env, jobject obj,
                                                                     jboolean narrowWindow) {
    Preferences::getInstance().setNarrowWindow(narrowWindow);
}

extern "C" JNIEXPORT void JNICALL
Java_fr_intuite_rtlsdrbridge_RtlSdrBridgeWrapper_nativeSetWwThresholdLVL2(JNIEnv *env, jobject obj,
                                                                        jfloat wwThresholdLVL2) {
    Preferences::getInstance().setWwThresholdLVL2(wwThresholdLVL2);
}

extern "C" JNIEXPORT void JNICALL
Java_fr_intuite_rtlsdrbridge_RtlSdrBridgeWrapper_nativeSetWwThresholdLVL3(JNIEnv *env, jobject obj,
                                                                        jfloat wwThresholdLVL3) {
    Preferences::getInstance().setWwThresholdLVL3(wwThresholdLVL3);
}

extern "C" JNIEXPORT jboolean JNICALL
Java_fr_intuite_rtlsdrbridge_RtlSdrBridgeWrapper_nativeInitParameters(
        JNIEnv *env,
        jobject /* this */,
        jlong centerFrequency,
        jlong sampleRate,
        jint samplesPerReading,
        jint freqFocusRangeKhz,
        jint gain,
        jlong refreshFFTMs,
        jlong refreshPeakMs,
        jlong refreshSignalStrengthMs,
        jboolean isMuted,
        jboolean dynamicThreshold,
        jboolean narrowWindow,
        jfloat wwThresholdLVL2,
        jfloat wwThresholdLVL3
) {
    Preferences &prefs = Preferences::getInstance();
// Initialize Preferences
    prefs.initialize(
            centerFrequency,
            sampleRate,
            samplesPerReading,
            freqFocusRangeKhz,
            gain,
            refreshFFTMs,
            refreshPeakMs,
            refreshSignalStrengthMs,
            isMuted,
            dynamicThreshold,
            narrowWindow,
            wwThresholdLVL2,
            wwThresholdLVL3
    );
    return true;
}

extern "C" JNIEXPORT jintArray JNICALL
Java_fr_intuite_rtlsdrbridge_RtlSdrBridgeWrapper_nativeGetTunerGains(JNIEnv *env, jobject obj) {
    if (dev == nullptr) {
        return nullptr;
    }

    int num_gains = rtlsdr_get_tuner_gains(dev, nullptr);
    if (num_gains <= 0) {
        return nullptr;
    }

    int *gains = new int[num_gains];
    rtlsdr_get_tuner_gains(dev, gains);

    jintArray result = env->NewIntArray(num_gains);
    env->SetIntArrayRegion(result, 0, num_gains, gains);

    delete[] gains;
    return result;
}