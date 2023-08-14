#ifdef _WIN32
// If you're running on Windows don't use MinGW stdin/stdout !!! These originate
// for C99/GNU99 compatibility. This application should be built with at least C11.
#define __USE_MINGW_ANSI_STDIO 0
#else
#define _WIN32 0
#endif

#include "fft.h"
#include "testFFT.h"

#include <math.h>
#include <complex.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

// Prints complex values to terminal
static bool PRINT_VALUES = true;

/**
 * Generates a pseudo random floating point number between 0 & given upper boundary
 * 
 * @param upperBoundary Upper boundary
 * @return float Floating point number between 0 & given upper boundary
 */
static inline float randomFloat(float upperBoundary) {
    return (float)rand()/(float)(RAND_MAX) * upperBoundary;
}

float compare(enum FFTAlgorithm algorithm, char algorithmName[], const float complex* testData, size_t size) {
    if(PRINT_VALUES) printf("Start comparison for %s algorithm...\n", algorithmName);
    clock_t start, stop;
    struct FFTResult *output, *outputInverse;
    if(PRINT_VALUES) printf("\nForward transformation (%s):\n\n", algorithmName);

    start = clock();
    output = fft(algorithm, testData, size);
    stop = clock();
    float duration = (float)(stop - start) / CLOCKS_PER_SEC;
    if(PRINT_VALUES) {
        for(size_t i=0; i<output->size; i++) {
            printf("TIME: (%f %+fj) ==> FREQUENCY: (%f %+fj)\n", crealf(testData[i]), cimagf(testData[i]), crealf(output->data[i]), cimagf(output->data[i]));
        }
        printf("\nInverse transformation (%s):\n\n", algorithmName);
    }

    start = clock();
    outputInverse = inverseFFT(algorithm, output->data, output->size);
    stop = clock();
    duration += (float)(stop - start) / CLOCKS_PER_SEC;
    if(PRINT_VALUES) {
        for(size_t i=0; i<outputInverse->size; i++) {
            printf("FREQUENCY: (%f %+fj) ==> TIME: (%f %+fj)\n", crealf(output->data[i]), cimagf(output->data[i]), crealf(outputInverse->data[i]), cimagf(outputInverse->data[i]));
        }
        printf("\nCompare original input with forward & backwards transformed output\n\n");

        for(size_t i=0; i<size; i++) {
            printf("ORIGINAL: (%f %+fj) ==> FORWARDS & BACKWARDS: (%f %+fj)\n", crealf(testData[i]), cimagf(testData[i]), crealf(outputInverse->data[i]), cimagf(outputInverse->data[i]));
            printf("DIFFERENCE: %f\n", cabsf(testData[i] - outputInverse->data[i]));
        }
    }

    deleteFFTResult(output);
    deleteFFTResult(outputInverse);
    return duration;
}

void forwardAndBackwardComparison(size_t noSamples, float upperBoundary) {
    // Initialize random generator
    srand((unsigned int)time(NULL));

    float complex* randomData = createDataArray(noSamples, noSamples);
    for (size_t i=0; i<noSamples; i++) {
        randomData[i] = randomFloat(upperBoundary) + randomFloat(upperBoundary) * I;
    }

    if(PRINT_VALUES) printf("\nStart now comparison for forward & backward transformation...\n");
    float time_dft = compare(BY_DEFINITION, "DFT", randomData, noSamples);
    if(PRINT_VALUES) printf("\n###############################################################################################\n\n");
    float time_cooley_tukey_recursive = compare(COOLEY_TUKEY_RECURSIVE, "Cooley-Tukey (recursive)", randomData, noSamples);
    if(PRINT_VALUES) printf("\n###############################################################################################\n\n");
    float time_cooley_tukey_iterative = compare(COOLEY_TUKEY_ITERATIVE, "Cooley-Tukey (iterative)", randomData, noSamples);
    if(PRINT_VALUES) printf("\n###############################################################################################\n\n");
    float time_bluestein = compare(BLUESTEIN, "Bluestein", randomData, noSamples);
    //if(PRINT_VALUES) printf("\n###############################################################################################\n\n");
    //float time_chirp_z = compare(CHIRP_Z, "Chirp-Z", randomData, size);

    printf("\n================================================================\n");
    printf("Duration (in sec.) for %zu samples\n", noSamples);
    printf("%-50s %.10fs\n", "Time for DFT:", time_dft);
    printf("%-50s %.10fs\n", "Time for Cooley-Tukey (recursive):", time_cooley_tukey_recursive);
    printf("%-50s %.10fs\n", "Time for Cooley-Tukey (iterative):", time_cooley_tukey_iterative);
    printf("%-50s %.10fs\n", "Time for Bluestein:", time_bluestein);
    //printf("%-50s %.10fs\n", "Time Chirp-Z:", time_chirp_z);
    printf("================================================================\n\n");
    free(randomData);
}

int main(int, char**) {

    // REMEMBER: To actually compare all (with & without power of 2) algorithms, we need arrays of size of power of 2!!!
    const size_t noSamples = 4096;

    float upperBoundary = 20.0f;

    // Setup buffer
    if(_WIN32 && PRINT_VALUES) {
        // Windows terminal/console is redicilous slow, therefore reduce printf calls by increasing the buffer
        char buffer[65536];
        setvbuf(stdout, buffer, _IOFBF, 65536);
    }

    if(testFFTAlgorithms()) {
        forwardAndBackwardComparison(noSamples, upperBoundary);
    }

    printf("Done...\n");
}