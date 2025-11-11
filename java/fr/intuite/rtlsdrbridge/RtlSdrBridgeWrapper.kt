package fr.intuite.rtlsdrbridge

object RtlSdrBridgeWrapper {

    // LogListener interface for handling native log messages
    fun interface LogListener {
        fun onLog(message: String)
    }

    private var logListener: LogListener? = null

    init {
        System.loadLibrary("rtl-sdr-bridge-android-lib")
    }

    // Set the LogListener
    fun setLogListener(listener: LogListener?) {
        logListener = listener
    }

    // Static method called from C++ to forward log messages
    @JvmStatic
    fun logFromNative(message: String) {
        logListener?.onLog(message) ?: android.util.Log.d("RtlSdrBridge", "C++: $message")
    }

    external fun nativeInitParameters(
        centerFrequency: Long,
        sampleRate: Long,
        samplesPerReading: Int,
        freqFocusRangeKhz: Int,
        gain: Int,
        refreshFFTMs: Long,
        refreshPeakMs: Long,
        refreshSignalStrengthMs: Long,
        ssbGain: Float,
        dynamicThreshold: Boolean,
        narrowWindow: Boolean,
        wwThresholdLVL2: Float,
        wwThresholdLVL3: Float
    ): Boolean

    external fun nativeInitRTL(
        fd: Int,
        path: String?
    ): Boolean

    external fun nativeReadAsync(
        fftCallback: (FloatArray) -> Unit,
        signalStrengthCallback: (Int) -> Unit,
        peakCallback: (Float) -> Unit,
        peakNormalizedCallback: (Float) -> Unit,
        peakFrequencyCallback: (Long) -> Unit,
        pcmCallback: (ShortArray) -> Unit
    )

    external fun nativeCancelAsync()

    external fun nativeCloseRTL()

    external fun nativeSetFrequency(frequency: Long)

    external fun nativeSetSampleRate(sampleRate: Long)

    external fun nativeSetGain(gain: Int)

    external fun nativeSetSamplesPerReading(samplesPerReading: Int)

    external fun nativeSetFrequencyFocusRange(frequencyFocusRange: Int)

    external fun nativeSetRefreshFFTMs(refreshFFTMs: Long)

    external fun nativeSetRefreshPeakMs(refreshPeakMs: Long)

    external fun nativeSetRefreshSignalStrengthMs(refreshSignalStrengthMs: Long)

    external fun nativeSetSsbGain(ssbGain: Float)

    external fun nativeSetNarrowWindow(narrowWindow: Boolean)

    external fun nativeSetDynamicThreshold(dynamicThreshold: Boolean)

    external fun nativeSetWwThresholdLVL2(wwThresholdLVL2: Float)

    external fun nativeSetWwThresholdLVL3(wwThresholdLVL3: Float)

    external fun nativeGetTunerGains(): IntArray?
}