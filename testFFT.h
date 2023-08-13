/**
 * API of some test cases for the FFT algorithms
 * 
 * @version 1.0
 * @date 2023-07-10
 * @file testFFT.h
 * @author Richard Saeuberlich (richard.saeuberlich@stud.tu-darmstadt.de)
 */
#pragma once

#include <stddef.h>
#include <stdbool.h>
#include <complex.h>

/**
 * Create an array which will be used to hold data in form of complex floating point values.
 * With this approach you may enter 
 * 
 * If 'realSize' does not equal 'size', the different entries will be zero padded.
 * 
 * Example: Size = 1000 will be rounded to 1024 = 2^10.
 * 
 * @param size Number of samples/data points
 * @param realSize Number of samples/data points with zero padding
 * @return float* Array with results
 */
float complex* createDataArray(size_t size, size_t realSize);

/**
 * Starts some tests for implemented FFT algorithms.
 */
bool testFFTAlgorithms();

/**
 * Makes a DFT & and FFT (Cooley-Tukey) with preset values and shows/prints,
 * the difference between these 2 algorithms.
 */
void showDifferenceDFTandCooleyTukeyRecursive();

/**
 * Makes a DFT & and FFT (Cooley-Tukey) with preset values and shows/prints,
 * the difference between these 2 algorithms.
 */
void showDifferenceDFTandCooleyTukeyIterative();

/**
 * Makes a 2 FFTs (Cooley-Tukey recursive & Cooley-Tukey iterative) with preset
 * values and shows/prints, the difference between these 2 algorithms.
 */
void showDifferenceCooleyTukeyRecursiveAndCooleyTukeyIterative();