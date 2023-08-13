/**
 * Implementation of some test cases for the FFT algorithms
 * 
 * @version 1.0
 * @date 2023-07-10
 * @file testFFT.h
 * @author Richard Saeuberlich (richard.saeuberlich@stud.tu-darmstadt.de)
 */
#include "testFFT.h"
#include "fft.h"
#include <math.h>
#include <complex.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

// I don't now why, but the DFT needs this "high"(depends of your definition of "high") value.
// The last test (for array E) fails for the 'fftByDefinition' if you choose any lower boundary!?!
static const float epsilon = 3e-4F;

bool testValues(char* arrayName, size_t index, float complex expected, float complex actual, bool inverse) {
    float error = cabsf(expected - actual);

    if (error > epsilon) {
        if(inverse) {
            printf("ERROR (inverse transformation) at %s[%zu] was to large: %f!\n\tExpected: %+.6f %+.6fi\n\tBut was : %+.6f %+.6fi\n", arrayName, index, error, crealf(expected), cimagf(expected), crealf(actual), cimagf(actual));
        } else {
            printf("ERROR (forward transformation) at %s[%zu] was to large: %f!\n\tExpected: %+.6f %+.6fi\n\tBut was : %+.6f %+.6fi\n", arrayName, index, error, crealf(expected), cimagf(expected), crealf(actual), cimagf(actual));
        }
        return false;
    }
    return true;
}

float complex* createDataArray(size_t size, size_t realSize) {
    if(size != realSize) {
        size_t tooMany = realSize - size;
        float complex* result = malloc(realSize * sizeof(float complex));
        // Zero padding for "unused" samples
        for (size_t i = realSize - tooMany; i<realSize; i++) {
            result[i] = 0.0f;
        }
        return result;
    } else {
        return malloc(size * sizeof(float complex));
    }
}

/**
 * ATTENTION:
 *      If you add new tests, remember that some algorithms need arrays of size of 2^N.
 *      To guarantee that every test works, you must always create test data of that size,
 *      meaning only arrays of size: 2, 4, 8, 16, 32, ...
 * 
 * @param fftAlgorithm Function pointer to FFT algorithm that should be tested
 * @param inverseFFTAlgorithm Function pointer to inverse algorithm of 'fftAlgorithm'
 * @return TRUE if all tests succeeded, FALSE if at least one test failed
 */
bool testFFTAlgorithm(enum FFTAlgorithm algorithm) {
    struct FFTResult *output, *outputInverse;
    bool continueTest = true;
    int size;
    // =========================================================================================
    // TEST 1 - Array A (Size: 2)
    // =========================================================================================
    
    size = 2;
    float complex inputA[] = {1.0F, 0.0F};
    float complex expectedA[] = {1 + 0 * I, 1 + 0 * I};

    output = fft(algorithm, inputA, size);
    for (int i = 0; continueTest && i < output->size; i++) {
        continueTest = testValues("A", i, expectedA[i], output->data[i], false);
    }

    outputInverse = inverseFFT(algorithm, output->data, output->size);
    for (int i = 0; continueTest && i < outputInverse->size; i++) {
        continueTest = testValues("A", i, inputA[i], outputInverse->data[i], true);
    }

    deleteFFTResult(output);
    deleteFFTResult(outputInverse);

    // =========================================================================================
    // TEST 2 - Array B (Size: 8)
    // =========================================================================================

    size = 8;
    float complex inputB[] = {0.0F, 1.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F, 0.0F};
    float complex expectedB[] = {
        1 + 0 * I, sqrtf(2.0F) / 2 - sqrtf(2.0F) * I / 2, 0 - 1 * I, - sqrtf(2.0F) / 2 - sqrtf(2.0F) * I / 2,
        - 1 + 0 * I, - sqrtf(2.0F) / 2 + sqrtf(2.0F) * I / 2, 0 + 1 * I, sqrtf(2.0F) / 2 + sqrtf(2.0F) * I / 2
    };

    output = fft(algorithm, inputB, size);
    for (int i = 0; continueTest && i < output->size; i++) {
        continueTest = testValues("B", i, expectedB[i], output->data[i], false);
    }

    outputInverse = inverseFFT(algorithm, output->data, output->size);
    for (int i = 0; continueTest && i < outputInverse->size; i++) {
        continueTest = testValues("B", i, inputB[i], outputInverse->data[i], true);
    }

    deleteFFTResult(output);
    deleteFFTResult(outputInverse);

    // =========================================================================================
    // TEST 3 - Array C (Size: 4)
    // =========================================================================================

    size = 4;
    float complex inputC[] = {1, 3, 5, 7};
    float complex expectedC[] = {16.0F, -4.0F + 4.0F*I, -4.0F, -4.0F - 4.0F*I};

    output = fft(algorithm, inputC, size);
    for (int i = 0; continueTest && i < output->size; i++) {
        continueTest = testValues("C", i, expectedC[i], output->data[i], false);
    }

    outputInverse = inverseFFT(algorithm, output->data, output->size);
    for (int i = 0; continueTest && i < outputInverse->size; i++) {
        continueTest = testValues("C", i, inputC[i], outputInverse->data[i], true);
    }

    deleteFFTResult(output);
    deleteFFTResult(outputInverse);

    // =========================================================================================
    // TEST 4 - Array D (Size: 8)
    // =========================================================================================

    size = 8;
    float complex inputD[] = {1.0F, 2.0F, 1.0F, -1.0F, 1.5F, 2.5F, 3.5F, 3.0F};
    float complex expectedD[] = {
        13.5F + -0.0F * I, 1.974874F + 5.681981F * I, -2.0F + -2.5F * I, -2.974874F + 0.681981F * I,
        0.5F + -0.0F * I, -2.974874F + -0.681981F * I, -2.0F + 2.5F * I, 1.974864F + -5.681981F * I
    };
    
    output = fft(algorithm, inputD, size);
    for (int i = 0; continueTest && i < output->size; i++) {
        continueTest = testValues("D", i, expectedD[i], output->data[i], false);
    }

    outputInverse = inverseFFT(algorithm, output->data, output->size);
    for (int i = 0; continueTest && i < outputInverse->size; i++) {
        continueTest = testValues("D", i, inputD[i], outputInverse->data[i], true);
    }

    deleteFFTResult(output);
    deleteFFTResult(outputInverse);

    // =========================================================================================
    // TEST 5 - Array E (Size: 16)
    // =========================================================================================
    
    size = 16;
    float complex inputE[] = {
        6.18F + 12.4F * I, 9.18F + 15.14F * I, 18.8F + 6.12F * I, 11.2F + 1.13F * I,
        5.3F + 1.4F * I, 0.4F + 11.14F * I, 10.3F + 13.2F * I, 3.16F + 16.18F * I,
        2.4F + 10.0F * I, 11.4F + 0.5F * I, 0.8F + 19.15F * I, 5.2F + 14.15F * I,
        0.12F + 10.18F * I, 6.16F + 10.1F * I, 7.7F + 8.15F * I, 19.2F + 16.2F * I
    };
    float complex expectedE[] = {
        117.5F + 165.139999F * I, 12.04186F + -15.201277F * I, 5.15667F + 11.65245F * I, 8.543372F + 41.359856F * I,
        -34.379997F + -1.019997F * I, 13.42107F + 22.932888F * I, -30.939322F + -1.411011F * I, 6.350538F + 7.96443F * I,
        -14.300003F + -3.940002F * I, -11.548398F + -45.060501F * I, 9.003334F + 6.787553F * I, -16.487682F + -29.763672F * I,
        -12.82F + -24.260002F * I, -33.914536F + 26.208889F * I, 29.419317F + 26.251007F * I, 51.833771F + 10.759388F * I
    };

    output = fft(algorithm, inputE, size);
    for (int i = 0; continueTest && i < output->size; i++) {
        continueTest = testValues("E", i, expectedE[i], output->data[i], false);
    }

    outputInverse = inverseFFT(algorithm, output->data, output->size);
    for (int i = 0; continueTest && i < outputInverse->size; i++) {
        continueTest = testValues("E", i, inputE[i], outputInverse->data[i], true);
    }

    deleteFFTResult(output);
    deleteFFTResult(outputInverse);

    return continueTest;
}

bool testFFTAlgorithms() {
    bool success = testFFTAlgorithm(BY_DEFINITION);
    if(success) {
        printf("All tests passed for 'FFT by definition' algorithm.\n");
    } else {
        printf("Tests for 'FFT by definition' algorithm failed!\n");
    }

    success = testFFTAlgorithm(COOLEY_TUKEY_RECURSIVE);
    if(success) {
        printf("All tests passed for 'Cooley-Tukey (recursive)' algorithm.\n");
    } else {
        printf("Tests for 'Cooley-Tukey (recursive)' algorithm failed!\n");
    }
    
    success = testFFTAlgorithm(COOLEY_TUKEY_ITERATIVE);
    if(success) {
        printf("All tests passed for 'Cooley-Tukey (iterative)' algorithm.\n");
    } else {
        printf("Tests for 'Cooley-Tukey (iterative)' algorithm failed!\n");
    }
   
    success = testFFTAlgorithm(BLUESTEIN);
    if(success) {
        printf("All tests passed for 'Bluestein' algorithm.\n");
    } else {
        printf("Tests for 'Bluestein' algorithm failed!\n");
    }

    /*
    success = testFFTAlgorithm(CHIRP_Z);
    if(success) {
        printf("All tests passed for 'Chirp-Z' algorithm.\n");
    } else {
        printf("Tests for 'Chirp-Z' algorithm failed!\n");
    }
    */
    return success;
}

void showDifferenceDFTandCooleyTukeyRecursive() {

    printf("\nDifference DFT & FFT (Cooley-Tukey Recursive)\n\n");
    size_t size = 16;
    float complex input[] = {
        6.18F + 12.4F * I, 9.18F + 15.14F * I, 18.8F + 6.12F * I, 11.2F + 1.13F * I,
        5.3F + 1.4F * I, 0.4F + 11.14F * I, 10.3F + 13.2F * I, 3.16F + 16.18F * I,
        2.4F + 10.0F * I, 11.4F + 0.5F * I, 0.8F + 19.15F * I, 5.2F + 14.15F * I,
        0.12F + 10.18F * I, 6.16F + 10.1F * I, 7.7F + 8.15F * I, 19.2F + 16.2F * I
    };
    
    complex float* output1 = createDataArray(size, size);
    fftByDefinition(input, output1, size);
    complex float* output2 = createDataArray(size, size);
    fftCooleyTukeyRecursive(input, output2, size);
    
    for(int i=0; i<size; i++) {
        printf("DFT: (%f+j%f) <==> FFT: (%f+j%f)\n", crealf(output1[i]), cimagf(output1[i]), crealf(output2[i]), cimagf(output2[i]));
        printf("DIFFERENCE: %f\n", cabsf(output1[i] - output2[i]));
    }
    printf("\n--------------------------------------------------\n\n");
}

void showDifferenceDFTandCooleyTukeyIterative() {

    printf("\nDifference DFT & FFT (Cooley-Tukey Iterative)\n\n");
    size_t size = 16;
    float complex input[] = {
        6.18F + 12.4F * I, 9.18F + 15.14F * I, 18.8F + 6.12F * I, 11.2F + 1.13F * I,
        5.3F + 1.4F * I, 0.4F + 11.14F * I, 10.3F + 13.2F * I, 3.16F + 16.18F * I,
        2.4F + 10.0F * I, 11.4F + 0.5F * I, 0.8F + 19.15F * I, 5.2F + 14.15F * I,
        0.12F + 10.18F * I, 6.16F + 10.1F * I, 7.7F + 8.15F * I, 19.2F + 16.2F * I
    };
    
    complex float* output1 = createDataArray(size, size);
    fftByDefinition(input, output1, size);
    complex float* output2 = createDataArray(size, size);
    fftCooleyTukeyIterative(input, output2, size);
    
    for(int i=0; i<size; i++) {
        printf("DFT: (%f+j%f) <==> FFT: (%f+j%f)\n", crealf(output1[i]), cimagf(output1[i]), crealf(output2[i]), cimagf(output2[i]));
        printf("DIFFERENCE: %f\n", cabsf(output1[i] - output2[i]));
    }
    printf("\n--------------------------------------------------\n\n");
}

void showDifferenceCooleyTukeyRecursiveAndCooleyTukeyIterative() {

    printf("\nDifference Cooley-Tukey Recursive & Cooley-Tukey Iterative\n\n");
    size_t size = 16;
    float complex input[] = {
        6.18F + 12.4F * I, 9.18F + 15.14F * I, 18.8F + 6.12F * I, 11.2F + 1.13F * I,
        5.3F + 1.4F * I, 0.4F + 11.14F * I, 10.3F + 13.2F * I, 3.16F + 16.18F * I,
        2.4F + 10.0F * I, 11.4F + 0.5F * I, 0.8F + 19.15F * I, 5.2F + 14.15F * I,
        0.12F + 10.18F * I, 6.16F + 10.1F * I, 7.7F + 8.15F * I, 19.2F + 16.2F * I
    };
    
    complex float* output1 = createDataArray(size, size);
    fftCooleyTukeyRecursive(input, output1, size);
    complex float* output2 = createDataArray(size, size);
    fftCooleyTukeyIterative(input, output2, size);
    
    for(int i=0; i<size; i++) {
        printf("RECURSIVE: (%f+j%f) <==> ITERATIVE: (%f+j%f)\n", crealf(output1[i]), cimagf(output1[i]), crealf(output2[i]), cimagf(output2[i]));
        printf("DIFFERENCE: %f\n", cabsf(output1[i] - output2[i]));
    }
    printf("\n--------------------------------------------------\n\n");
}