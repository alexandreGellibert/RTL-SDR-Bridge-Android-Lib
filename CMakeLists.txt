cmake_minimum_required(VERSION 3.22.1)

project("rtl-sdr-bridge-android-lib")

# Force 16 KB segment alignment
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-z,max-page-size=0x4000")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,-z,max-page-size=0x4000")

# Build libusb as a shared library
add_subdirectory(libusb)

# fftw3f (static)
add_library(fftw3f STATIC IMPORTED)
set_target_properties(fftw3f PROPERTIES
        IMPORTED_LOCATION "${CMAKE_SOURCE_DIR}/libs/${ANDROID_ABI}/libfftw3f.a"
        INTERFACE_INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/include/fftw3"
)

# App shared library
add_library(${CMAKE_PROJECT_NAME} SHARED
        # List C/C++ source files with relative paths to this CMakeLists.txt.
        src/rtl-sdr-bridge-java.cpp
        # src/ssb_demod.cpp
        src/ssb_demod_high.cpp
        # src/ssb_demod_boost.cpp
        rtl-sdr/src/librtlsdr.c
        rtl-sdr/src/tuner_e4k.c
        rtl-sdr/src/tuner_fc0012.c
        rtl-sdr/src/tuner_fc0013.c
        rtl-sdr/src/tuner_fc2580.c
        rtl-sdr/src/tuner_r82xx.c
        rtl-sdr/src/convenience/convenience.c)

target_include_directories(${CMAKE_PROJECT_NAME} PRIVATE
        ${CMAKE_SOURCE_DIR}/rtl-sdr/include
        ${CMAKE_SOURCE_DIR}/libusb/libusb
        ${CMAKE_SOURCE_DIR}/libusb/android
        ${CMAKE_SOURCE_DIR}
)

target_compile_options(${CMAKE_PROJECT_NAME} PRIVATE -w)

target_link_libraries(${CMAKE_PROJECT_NAME}
        PRIVATE
        android
        log
        fftw3f
        libusb
)