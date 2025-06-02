# RTL-SDR Bridge for Java/Android lib
A C++ library for RTL-SDR and FFTW with JNI bindings for Android.

## Build Instructions
1. Install CMake and Android NDK.
2. Run: `cmake . && make`
3. Copy the generated .so files to your Android app's jniLibs.

## Usage
```java
import com.intuite.rtlsdrbridge.RtlSdrWrapper;
RtlSdrWrapper wrapper = new RtlSdrWrapper();
wrapper.init();
wrapper.readAsync(callbacks);
