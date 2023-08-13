/**
 * API of Fast Fourier Transformation (FFT) in 5 flavours:
 *      1.) FFT by definition algorithm (actually DFT)
 *      2.) Cooley–Tukey FFT algorithm (recursive)
 *      3.) Cooley–Tukey FFT algorithm (iterative)
 *      4.) Bluestein FFT algorithm
 *      5.) Chirp-Z FFT algorithm (!!!Does currently NOT work!!!)
 * 
 * @version 1.0
 * @date 2023-07-10
 * @file fft.h
 * @author Richard Saeuberlich (richard.saeuberlich@stud.tu-darmstadt.de)
 */
#pragma once

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>
#include <complex.h>

/**
 * Holds actual array (with data) & size of this array
 */
struct FFTResult {
    float complex* data;
    size_t size;
};

/**
 * Bit significance:
 *      - MSB (Most significant bit): Most left one
 *      - LSB (Least significant bit): Most right one
 */
enum BitSignificance {
    MSB, LSB
};

/**
 * Enumeration of all implemented FFT algorithms
 */
enum FFTAlgorithm {
    BY_DEFINITION,
    COOLEY_TUKEY_RECURSIVE,
    COOLEY_TUKEY_ITERATIVE,
    BLUESTEIN,
    CHIRP_Z
};

static const bool CHECK_ARGUMENTS = true;
static const bool PRINT_REVERSE_INDEXES = false;
static const enum BitSignificance SIGNIFICANCE = MSB;

/**
 * This is not really a FFT but the implement of the
 * Discrete Fourier Transformation (DFT) by definition.
 * 
 * @param dataIN Samples/Data array of size 'size'
 * @param dataOUT Result array of size 'size'
 * @param size Number of samples/data points
 * @see https://en.wikipedia.org/wiki/Fast_Fourier_transform#Definition
 */
void fftByDefinition(const float complex* dataIN, float complex* dataOUT, size_t size);

/**
 * This is not really a FFT but the implement of the
 * Discrete Fourier Transformation (DFT) by definition.
 * 
 * @param dataIN Samples/Data array of size 'size'
 * @param dataOUT Result array of size 'size'
 * @param size Number of samples/data points
 * @see https://en.wikipedia.org/wiki/Fast_Fourier_transform#Definition
 */
void inverseFFTByDefinition(const float complex* dataIN, float complex* dataOUT, size_t size);

/**
 * This is the classical FFT implementation after the Cooley–Tukey algorithm.
 * 
 * @param dataIN Samples/Data array of size 'size' (must be of length of power of 2)
 * @param dataOUT Result array of size 'size' (must be of length of power of 2)
 * @param size Number of samples/data points (must be power of 2)
 * @see https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
 */
void fftCooleyTukeyRecursive(const float complex* dataIN, float complex* dataOUT, size_t size);

/**
 * This is the inverse FFT of the 'fftCooleyTukeyRecursive'.
 * 
 * @param dataIN Samples/Data array of size 'size' (must be of length of power of 2)
 * @param dataOUT Result array of size 'size' (must be of length of power of 2)
 * @param size Number of samples/data points (must be power of 2)
 * @see https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
 */
void inverseFFTCooleyTukeyRecursive(const float complex* dataIN, float complex* dataOUT, size_t size);

/**
 * This Cooley-Tukey FFT implementation works NOT recursive, but with a reordering
 * of the resulting array. In this specific case the radix-2 decimation-in-time (DIT) is used.
 * 
 * @param dataIN Samples/Data array of size 'size' (must be of length of power of 2)
 * @param dataOUT Result array of size 'size' (must be of length of power of 2)
 * @param size Number of samples/data points (must be power of 2)
 * @see https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
 */
void fftCooleyTukeyIterative(const float complex* dataIN, float complex* dataOUT, size_t size);

/**
 * This is the inverse FFT of the 'fftCooleyTukeyIterative'.
 * 
 * @param dataIN Samples/Data array of size 'size' (must be of length of power of 2)
 * @param dataOUT Result array of size 'size' (must be of length of power of 2)
 * @param size Number of samples/data points (must be power of 2)
 * @see https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
 */
void inverseFFTCooleyTukeyIterative(const float complex* dataIN, float complex* dataOUT, size_t size);

/**
 * The Bluestein FFT algorithm, needs some zero initialized arrays.
 * This is a quick/efficient way to fill a given array with zeros.
 * 
 * @param data Array that will hold only zeros afterwards of size 'size'
 * @param size Size of given array
 */
void resetArray(float complex* data, size_t size);

/**
 * Main advantage of Bluestein algorithm is, you don't need an array of size of power 2
 * 
 * @param dataIN Samples/Data array of size 'size'
 * @param dataOUT Result array of size 'size'
 * @param twiddleFactors Array of size 'size', array doesn't need to be initialized
 * @param a Array of size 'powerOfTwoSize', all entries must be initialized with 0
 * @param b Array of size 'powerOfTwoSize', all entries must be initialized with 0
 * @param c Array of size 'powerOfTwoSize', array doesn't need to be initialized with zeros
 * @param d Array of size 'powerOfTwoSize', array doesn't need to be initialized with zeros
 * @param powerOfTwoSize Next power-of-2 such that following equation is true: 'powerOfTwoSize' >= 'size' * 2 + 1
 * @param size Number of samples/data points
 * @see https://en.wikipedia.org/wiki/Chirp_Z-transform#Bluestein's_algorithm
 */
void fftBluestein(const float complex* dataIN, float complex* dataOUT, float complex* twiddleFactors, float complex* a, float complex* b, float complex* c, float complex* d, size_t powerOfTwoSize, size_t size);

/**
 * This is the inverse FFT of the 'fftBluestein'.
 * 
 * @param dataIN Samples/Data array of size 'size'
 * @param dataOUT Result array of size 'size'
 * @param twiddleFactors Array of size 'size', array doesn't need to be initialized for twiddle calculation
 * @param a Array of size 'powerOfTwoSize', all entries must be initialized with 0
 * @param b Array of size 'powerOfTwoSize', all entries must be initialized with 0
 * @param c Array of size 'powerOfTwoSize', array doesn't need to be initialized with zeros
 * @param d Array of size 'powerOfTwoSize', array doesn't need to be initialized with zeros
 * @param powerOfTwoSize Next power-of-2 such that following equation is true: 'powerOfTwoSize' >= 'size' * 2 + 1
 * @param size Number of samples/data points
 * @see https://en.wikipedia.org/wiki/Chirp_Z-transform#Bluestein's_algorithm
 */
void inverseFFTBluestein(const float complex* dataIN, float complex* dataOUT, float complex* twiddleFactors, float complex* a, float complex* b, float complex* c, float complex* d, size_t powerOfTwoSize, size_t size);

/**
 * Main advantage of Chirp-Z algorithm is, you don't need an array of size of power 2
 * 
 * @param dataIN Samples/Data array of size 'size'
 * @param dataOUT Result
 * @param chirps Empty array for chirp calculation
 * @param weights Empty array for weight calculation
 * @param size Number of samples/data points
 * @see https://en.wikipedia.org/wiki/Chirp_Z-transform
 */
void fftChirpZ(const float complex* dataIN, float complex* dataOUT, float complex* chirps, float complex* weights, size_t size);

/**
 * This is the inverse FFT of the 'fftChirpZ'
 * 
 * @param dataIN Samples/Data array of size 'size'
 * @param dataOUT Result array of size 'size'
 * @param chirps Array of size 'size' for chirp calculation
 * @param weights Array of size 'size' for weight calculation
 * @param size Number of samples/data points
 * @see https://en.wikipedia.org/wiki/Chirp_Z-transform
 */
void inverseFFTChirpZ(const float complex* dataIN, float complex* dataOUT, float complex* chirps, float complex* weights, size_t size);

// ==========================================================================
// Convenient API functions
// ==========================================================================

/**
 * Deletes given FFTResult struct
 * 
 * @param result FFTResult struct to be deleted
 */
void deleteFFTResult(struct FFTResult* result);

/**
 * Executes FFT with desired algorithm
 * 
 * NOTE: Use this convenient function call only if you're allowed to place more complex data structures on the stack without any problems.
 *       If you are not allowed to do this, you should call the individual FFT functions directly, with pre-initialized arrays.
 * 
 * @param algorithm Algorithm of choice
 * @param data Samples/Data array of size 'size'
 * @param size Number of samples/data points
 * @return Struct with new allocated array (resulting data) & size, or NULL if chosen FFT algorithm was unknown
 */
struct FFTResult* fft(enum FFTAlgorithm algorithm, const float complex* data, size_t size);

/**
 * Executes inverse FFT with desired algorithm
 * 
 * NOTE: Use this convenient function call only if you're allowed to place more complex data structures on the stack without any problems.
 *       If you are not allowed to do this, you should call the individual FFT functions directly, with pre-initialized arrays.
 * 
 * @param algorithm Algorithm of choice
 * @param data Samples/Data array of size 'size'
 * @param size Number of samples/data points
 * @return Struct with new allocated array (resulting data) & size, or NULL if chosen FFT algorithm was unknown
 */
struct FFTResult* inverseFFT(enum FFTAlgorithm algorithm, const float complex* data, size_t size);