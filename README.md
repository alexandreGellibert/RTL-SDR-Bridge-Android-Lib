# RTL-SDR Bridge for Android lib
A C++ library for RTL-SDR and FFTW with JNI bindings for Android.

The Wrapper is written in Kotlin but can be adapted for Java without too much trouble.

Our algorithm transforms the signal (FFT), detects peaks and signal quality too.

## How to install
1. Create an Android Project (Kotlin - Java 17) (if you use another Java version, you need to update the buid.gradle.ktx file)
2. Add the RTL-SDR-Bridge library as a git submodule:
```bash
git add .gitmodules rtl-sdr-lib
git commit -m "Add rtl-sdr-lib as a submodule"
git push origin main
```

3. Add RTL-SDR-Bridge as a dependance in your Android Project

app/build.gradle.kts
```
implementation(project(":RTL-SDR-Bridge-Android-Lib"))Add commentMore actions
project(":RTL-SDR-Bridge-Android-Lib")
```
settings.gradle.kts
```
include(":app", ":RTL-SDR-Bridge-Android-Lib")
```
gradle/libs.versions.toml
```
[versions]
kotlin = "1.9.22" #or another version

[libraries]
androidx-core-ktx = { group = "androidx.core", name = "core-ktx", version.ref = "coreKtx" }

[plugins]
kotlin-android = { id = "org.jetbrains.kotlin.android", version = "kotlin" }
```



## How to use

The Wrapper has a lot of methods to init and set the different parameters of the RTL-SDR : sample rate, gain, center frequency etc.
Once it starts reading, it will callback the Kotlin/Java through 4 callbacks:
1. fftCallback: (FloatArray) -> Unit = {} : A Float Array of the FFT of the sample processed
2. signalStrengthCallback: (Int) -> Unit  = {} : An indication if signal was found or not (0 no signal to 3 good signal)
3. peakCallback: (Float) -> Unit = {} : The strength of the signal at its Peak
4. peakFrequencyCallback: (Long) -> Unit = {} : The frequency of the signal at its Peak

```java
import fr.intuite.rtlsdrbridge.RtlSdrBridgeWrapper

private val rtlBridgeWrapper = RtlSdrBridgeWrapper

rtlBridgeWrapper.setLogListener { message -> // how you handle C++ logs }

success = rtlBridgeWrapper.nativeInitRTL(deviceFileDescriptor, devicePath)
if (success) {
  rtlBridgeWrapper.nativeReadAsync(fftCallback, signalStrengthCallback, peakCallback, peakFrequencyCallback)
}
```

## Open-source Libraries used
Thanks for all the great open-source projects that were used in this project:

FFTW3 https://github.com/FFTW/fftw3

RTL-SDR (modified) https://github.com/osmocom/rtl-sdr.git

libusb  https://github.com/libusb/libusb.git
