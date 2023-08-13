/**
 * Implementation of Fast Fourier Transformation (FFT) in 5 flavours:
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
#include "constants.h"
#include "fft.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

// ==========================================================================
// Utility functions
// ==========================================================================

/**
 * Checks if given arrays are NULL
 * 
 * @param dataIN Data input array to be tested
 * @param dataOUT Data output array to be tested
 * @return TRUE if at least one array is NULL, FALSE otherwise
 */
static bool dataArraysAreInvalid(const float complex* dataIN, const float complex* dataOUT) {
    return dataIN == NULL || dataOUT == NULL;
}

/**
 * Calculates the exponent of the nth root (N) of unity: e ^ +j * 2 * pi / N
 * 
 * In the context of a FFT, this is often called: twiddle factor.
 * 
 * @param size Amount of samples/data points
 * @param withSize If TRUE, amount of samples/data points will be considered. If FALSE, amount of samples/data points will be ignored.
 * @return float complex Exponent that may be used to calculate twiddles
 * @see https://en.wikipedia.org/wiki/Root_of_unity
 * @see https://en.wikipedia.org/wiki/Twiddle_factor
 */
static inline float complex calcForward(size_t size, bool withSize) {
    if(withSize) {
        return -I * 2.0f * PI / size;
    } else {
        return -I * 2.0f * PI;
    }
}

/**
 * Calculates the complex conjugated exponent of the nth root (N) of unity: e ^ -j * 2 * pi / N
 * 
 * In the context of a FFT, this is often called: twiddle factor.
 * 
 * @param size Amount of samples/data points
 * @param withSize If TRUE, amount of samples/data points will be considered. If FALSE, amount of samples/data points will be ignored.
 * @return float complex Exponent that may be used to calculate twiddles
 * @see https://en.wikipedia.org/wiki/Root_of_unity
 * @see https://en.wikipedia.org/wiki/Twiddle_factor
 */
static inline float complex calcInverse(size_t size, bool withSize) {
    return -1.0f * calcForward(size, withSize);
}

/**
 * Returns highest bit position for given 32 bit integer with desired bit significance.
 * Example for 37 = 0010 0101:
 *      - MSB: Highest bit = 2, because first 1 from left to right (->) is at position 2
 *      - LSB: Highest bit = 0, because first 1 from right to left (<-) is at position 0
 * 
 * @param size Number
 * @param significance Search for first bit from left (MSB) or right (LSB)
 * @return uint8_t Highest bit position
 */
static uint8_t findHighestBit(size_t size, enum BitSignificance significance) {
    if(significance == LSB) {
        for(uint8_t i=0; i<32; i++) {
            if((size >> i) & 1) {
                return i;
            }
        }
    } else {
        for(uint8_t i=31; i>=0; i--) {
            if((size >> i) & 1) {
                return i;
            }
        }
    }
}

/**
 * Some FFT algorithms need a sample size of power of 2. This function checks if given
 * size is of power of 2. If so, same size will be returned. If not, the size will be
 * rounded up to the next power of 2.
 * 
 * @param size Actual sample/data point size
 * @param significance MSB or LSB
 * @return size_t New size rounded to the next power of 2
 */
static size_t adjustSize(size_t size) {
    uint8_t highestBit = findHighestBit(size, SIGNIFICANCE);
    if(size > (1 << highestBit)) {
        // Round to next power of 2
        return 1 << (highestBit + 1);
    } else {
        // Exact fit
        return 1 << highestBit;
    }
}

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
static float complex* createDataArray(size_t size, size_t realSize) {
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
 * Squares given integer
 * 
 * @param number Number to be squared
 * @return uint32_t Result
 */
static inline size_t square(size_t number) {
    return number * number;
}

inline void resetArray(float complex* data, size_t size) {
    memset(data, 0, size*sizeof(float complex));
}

// ==========================================================================
// FFT by definition (a.k.a. DFT)
// @see https://en.wikipedia.org/wiki/Discrete_Fourier_transform
// ==========================================================================

void fftByDefinition(const float complex* dataIN, float complex* dataOUT, size_t size) {
    if(CHECK_ARGUMENTS) {
        if(dataArraysAreInvalid(dataIN, dataOUT)) {
            printf("Given arrays mustn't be NULL!!!");
            return;
        }
    }
    float complex coefficient = calcForward(size, true);
    for (size_t k = 0; k < size; k++) {
        dataOUT[k] = 0;
        for (size_t n = 0; n < size; n++) {
            dataOUT[k] += dataIN[n] * cexpf(coefficient * n * k);
        }
    }
}

void inverseFFTByDefinition(const float complex* dataIN, float complex* dataOUT, size_t size) {
    if(CHECK_ARGUMENTS) {
        if(dataArraysAreInvalid(dataIN, dataOUT)) {
            printf("Given arrays mustn't be NULL!!!");
            return;
        }
    }
    float complex coefficient = calcInverse(size, true);
    for (size_t k = 0; k < size; k++) {
        dataOUT[k] = 0;
        for (size_t n = 0; n < size; n++) {
            dataOUT[k] += dataIN[n] * cexpf(coefficient * n * k);
        }
        dataOUT[k] /= size;
    }
}

// ==========================================================================
// Cooley-Tukey recursive
// @see https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
// ==========================================================================

/**
 * Cooley-Tukey FFT with recursive algorithm & out-of-place array for the results
 * 
 * Remember:
 *      The pointer 'dataIN' is NOT at the beginning of the actual array, because of some pointer magic.
 *      This means dataIN[0] is NOT (except one time) the first value!!! Same goes for 'dataOUT' pointer.
 * 
 * @param dataIN 
 * @param dataOUT 
 * @param size 
 * @param offset 
 */
static void _fftCooleyTukeyRecursive(const float complex* dataIN, float complex* dataOUT, size_t size, size_t offset, const float complex coefficient) {
    if (size < 2) {
        // If size is 1, just pass through
        dataOUT[0] = dataIN[0];
    } else {
        size_t halvedSize = size >> 1;        // REMEMBER: x/2 -> x >> 1
        size_t doubledOffset = offset << 1;   // REMEMBER: x*2 -> x << 1
        
        // Mimic array splitting with some pointer magic...
        _fftCooleyTukeyRecursive(dataIN, dataOUT, halvedSize, doubledOffset, coefficient);
        _fftCooleyTukeyRecursive(dataIN + offset, dataOUT + halvedSize, halvedSize, doubledOffset, coefficient);
        
        // Reattach the 2 arrays into one
        register size_t newIndex;
        register float complex tmp1, tmp2;
        for (size_t k = 0; k < halvedSize; k++) {
            newIndex = k + halvedSize;
            tmp1 = dataOUT[k];
            tmp2 = cexpf(coefficient * k / size) * dataOUT[newIndex];

            dataOUT[k] = tmp1 + tmp2;
            dataOUT[newIndex] = tmp1 - tmp2;
        }
    }
}

void fftCooleyTukeyRecursive(const float complex* dataIN, float complex* dataOUT, size_t size) {
    if(CHECK_ARGUMENTS) {
        if(dataArraysAreInvalid(dataIN, dataOUT)) {
            printf("Given arrays mustn't be NULL!!!");
            return;
        }
        if(size != adjustSize(size)) {
            printf("Given size (%lu) is NOT a power of 2", size);
            return;
        }
    }
    _fftCooleyTukeyRecursive(dataIN, dataOUT, size, 1, calcForward(size, false));
}

void inverseFFTCooleyTukeyRecursive(const float complex* dataIN, float complex* dataOUT, size_t size) {
    if(CHECK_ARGUMENTS) {
        if(dataArraysAreInvalid(dataIN, dataOUT)) {
            printf("Given arrays mustn't be NULL!!!");
            return;
        }
        if(size != adjustSize(size)) {
            printf("Given size (%lu) is NOT a power of 2", size);
            return;
        }
    }
    _fftCooleyTukeyRecursive(dataIN, dataOUT, size, 1, calcInverse(size, false));
    for (size_t i = 0; i < size; ++i) {
        dataOUT[i] /= size;
    }
}

// ==========================================================================
// Cooley-Tukey iterative
// @see https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
// ==========================================================================

/**
 * Copies input data in bit-reversed order into output array.
 * This is necessary for the Cooley-Tukey iteration FFT algorithm.
 * 
 * @param dataIN Array with original ordered data
 * @param dataOUT Array that will hold the swapped data afterwards
 * @param size Size of given arrays
 * @see https://en.wikipedia.org/wiki/Butterfly_diagram#Radix-2_butterfly_diagram
 */
static void bitReverseCopy(const float complex* dataIN, float complex* dataOUT, const size_t size) {
    size_t forwards, backwards, groupSize;

    // First and last element stay the same!
    dataOUT[0] = dataIN[0];
    dataOUT[size - 1] = dataIN[size - 1];

    // Therefore start at 1 & stop at size-1!
    for (forwards = 1, backwards = size / 2; forwards < size - 1; forwards++) {
        // Swap only as long as you are in lower part (meaning first half) of array
        if (forwards <= backwards) {
            dataOUT[forwards] = dataIN[backwards];
            dataOUT[backwards] = dataIN[forwards];
        }
        groupSize = size / 2;
        /*
         * This while loop effectively performs a bit-reversal operation on "backwards".
         * The variable "groupSize" represents half an array. This halfs every iteration.
         * By repeatedly subtracting "groupSize" from "backwards", the loop traverses the bits of "backwards" in reverse order.
         */
        while (groupSize <= backwards) {
            backwards -= groupSize;
            groupSize /= 2;
        }
        // Goto next bit reversal index
        backwards += groupSize;
    }
}

/**
 * Cooley-Tukey FFT with iterative algorithm
 * 
 * @param dataIN Array with data
 * @param dataOUT Array that will hold the transformed data afterwards
 * @param size Size of given arrays
 * @param coefficient Precalculated unity root coefficiant
 * @see https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm#Data_reordering,_bit_reversal,_and_in-place_algorithms
 */
static void _fftCooleyTukeyIterative(const float complex* dataIN, float complex* dataOUT, const size_t size, const float complex coefficient) {
    bitReverseCopy(dataIN, dataOUT, size);

    size_t groupSize, currentGroupSize, i, j;
    // Twiddles are unity roots (@see https://en.wikipedia.org/wiki/Twiddle_factor)
    float complex twiddle, twiddleProduct, tmp;
    for (groupSize = 2; groupSize <= size; groupSize = groupSize << 1) { // REMEMBER: x*2 = x << 1
        twiddle = cexpf(coefficient / groupSize);
        for (i = 0; i < size; i += groupSize) {
            twiddleProduct = 1; // REMEMBER: e^0 = 1
            currentGroupSize = groupSize >> 1; // REMEMBER: x/2 = x >> 1
            for (j = 0; j < currentGroupSize; ++j) {
                tmp = dataOUT[i + j + currentGroupSize] * twiddleProduct;
                dataOUT[i + j + currentGroupSize] = dataOUT[i + j] - tmp;
                dataOUT[i + j] += tmp;
                twiddleProduct *= twiddle;
            }
        }
    }
}

void fftCooleyTukeyIterative(const float complex* dataIN, float complex* dataOUT, size_t size) {
    if(CHECK_ARGUMENTS) {
        if(dataArraysAreInvalid(dataIN, dataOUT)) {
            printf("Given arrays mustn't be NULL!!!");
            return;
        }
        if(size != adjustSize(size)) {
            printf("Given size (%lu) is NOT a power of 2", size);
            return;
        }
    }
    _fftCooleyTukeyIterative(dataIN, dataOUT, size, calcForward(size, false));
}

void inverseFFTCooleyTukeyIterative(const float complex* dataIN, float complex* dataOUT, size_t size) {
    if(CHECK_ARGUMENTS) {
        if(dataArraysAreInvalid(dataIN, dataOUT)) {
            printf("Given arrays mustn't be NULL!!!");
            return;
        }
        if(size != adjustSize(size)) {
            printf("Given size (%lu) is NOT a power of 2", size);
            return;
        }
    }
    _fftCooleyTukeyIterative(dataIN, dataOUT, size, calcInverse(size, false));
    for (size_t i = 0; i < size; ++i) {
        dataOUT[i] /= size;
    }
}

// ==========================================================================
// Bluestein FFT
// @see https://en.wikipedia.org/wiki/Chirp_Z-transform#Bluestein's_algorithm
// ==========================================================================

/**
 * Convolution in time domain is equal to a multiplication in frequency domain:
 *      dataOUT = (f * g) = IDFT(DFT(f) * DFT(g))
 * 
 * @param f Samples/Data points (of size 'powerOfTwoSize') of first signal
 * @param g Samples/Data points (of size 'powerOfTwoSize') of second signal
 * @param h Just a "working space" array
 * @param dataOUT Convolution of first & second signal
 * @param powerOfTwoSize Number of samples/data points (must be power of 2)
 * @see https://en.wikipedia.org/wiki/Convolution_theorem
 */
static void convolution(const float complex* f, const float complex* g, float complex* h, float complex* dataOUT, const size_t powerOfTwoSize) {
    fftCooleyTukeyIterative(f, dataOUT, powerOfTwoSize);
    fftCooleyTukeyIterative(g, h, powerOfTwoSize);

	for (size_t i = 0; i < powerOfTwoSize; i++) {
		h[i] *= dataOUT[i];
	}
    
    inverseFFTCooleyTukeyIterative(h, dataOUT, powerOfTwoSize);
}

static void _fftBluestein(const float complex* dataIN, float complex* dataOUT, float complex* twiddleFactors, float complex* a, float complex* b, float complex* c, float complex* d, const size_t powerOfTwoSize, const size_t size, bool inverse) {
    size_t temp;
    float complex twiddle;
    const float coefficient = (inverse ? 1.0f : -1.0f) * PI / size;
	for (size_t i = 0; i < size; i++) {
        temp = square(i) % ((size_t)size * 2u);
        twiddle = cexpf(coefficient * temp * I);
		twiddleFactors[i] = twiddle;
        a[i] = dataIN[i] * twiddle;
        b[i] = conjf(twiddle);
	}
	
	b[0] = twiddleFactors[0];
	for (size_t i = 1; i < size; i++) {
        b[powerOfTwoSize - i] = conjf(twiddleFactors[i]);
    }
    
	convolution(a, b, c, d, powerOfTwoSize);

    float tempSize = (inverse ? size : 1);
    for (size_t i = 0; i < size; i++) {
        dataOUT[i] = d[i] * twiddleFactors[i] / tempSize;
    }
}

void fftBluestein(const float complex* dataIN, float complex* dataOUT, float complex* twiddleFactors, float complex* a, float complex* b, float complex* c, float complex* d, size_t powerOfTwoSize, size_t size) {
    if(CHECK_ARGUMENTS) {
        if(dataArraysAreInvalid(dataIN, dataOUT)) {
            printf("Given arrays mustn't be NULL!!!");
            return;
        }
        if(powerOfTwoSize != adjustSize(powerOfTwoSize)) {
            printf("Given powerOfTwoSize (%lu) is NOT a power of 2", powerOfTwoSize);
            return;
        }
    }
    _fftBluestein(dataIN, dataOUT, twiddleFactors, a, b, c, d, powerOfTwoSize, size, false);
}

void inverseFFTBluestein(const float complex* dataIN, float complex* dataOUT, float complex* twiddleFactors, float complex* a, float complex* b, float complex* c, float complex* d, size_t powerOfTwoSize, size_t size) {
    if(CHECK_ARGUMENTS) {
        if(dataArraysAreInvalid(dataIN, dataOUT)) {
            printf("Given arrays mustn't be NULL!!!");
            return;
        }
        if(powerOfTwoSize != adjustSize(powerOfTwoSize)) {
            printf("Given powerOfTwoSize (%lu) is NOT a power of 2", powerOfTwoSize);
            return;
        }
    }
    _fftBluestein(dataIN, dataOUT, twiddleFactors, a, b, c, d, powerOfTwoSize, size, true);
}

// ==========================================================================
// Chirp-Z FFT (!!!Does currently NOT work!!!)
// @see https://en.wikipedia.org/wiki/Chirp_Z-transform
// ==========================================================================

/**
 * Calculates modulo with respect to given size
 * 
 * @param a 
 * @param size 
 * @return size_t 
 */
static size_t modulo(size_t a, size_t size) {
    size_t result = a % size;
    if (result < 0) {
        result += size;
    }
    return result;
}

/**
 * Chirps are twiddle factors. For some reason, people called it chirps instead of twiddles...?!?
 * 
 * Calculates sum from n=0 to n=size-1 of:
 *      e ^ (-j * 2pi * n^2 * omega / size)
 * 
 * @param chirps Array which holds the calculated chirps afterwards
 * @param omega ???
 * @param size Amount of samples
 * @see https://en.wikipedia.org/wiki/Chirp_Z-transform
 * @see https://en.wikipedia.org/wiki/Twiddle_factor
 */
static void calculateChirps(float complex* chirps, float omega, size_t size) {
    const float complex coefficient = calcForward(size, true) * omega;
    for (size_t n = 0; n < size; n++) {
        chirps[n] = cexpf(coefficient * square(n));
    }
}

/**
 * Calculates sum from k=0 to k=size-1 of:
 *      e ^ (-j * 2pi * k / size)
 * 
 * @param weight Array which holds the calculated weights afterwards
 * @param size Amount of samples
 * @see https://en.wikipedia.org/wiki/Chirp_Z-transform
 */
static void calculateWeights(float complex* weights, size_t size) {
    const float complex coefficient = calcForward(size, true);
    for (size_t k = 0; k < size; k++) {
        weights[k] = cexpf(coefficient * k);
    }
}

static void _fftChirpZ(const float complex* dataIN, float complex* dataOUT, float complex* chirps, float complex* weights, float omega, size_t size) {

    calculateChirps(chirps, omega, size);

    calculateWeights(weights, size);
    
    // Chirp-Z transformation
    size_t index;
    for (size_t k = 0; k < size; k++) {
        dataOUT[k] = 0;
        for (size_t n = 0; n < size; n++) {
            index = modulo(k * n, size);
            dataOUT[k] += dataIN[n] * chirps[index] * weights[modulo(-k * n, size)];
        }
    }
}

void fftChirpZ(const float complex* dataIN, float complex* dataOUT, float complex* chirps, float complex* weights, size_t size) {
    if(CHECK_ARGUMENTS) {
        if(dataArraysAreInvalid(dataIN, dataOUT)) {
            printf("Given arrays mustn't be NULL!!!");
            return;
        }
    }
    // TODO
    float omega = 0.1f; // What is omega?!?
    _fftChirpZ(dataIN, dataOUT, chirps, weights, omega, size);
}

void inverseFFTChirpZ(const float complex* dataIN, float complex* dataOUT, float complex* chirps, float complex* weights, size_t size) {
    if(CHECK_ARGUMENTS) {
        if(dataArraysAreInvalid(dataIN, dataOUT)) {
            printf("Given arrays mustn't be NULL!!!");
            return;
        }
    }
    // TODO
}

// ==========================================================================
// Convenient API functions
// ==========================================================================

void deleteFFTResult(struct FFTResult* result) {
    free(result->data);
    free(result);
}

struct FFTResult* fft(enum FFTAlgorithm algorithm, const float complex* data, size_t size) {
    size_t resultSize;
    float complex* dataOUT;
    switch (algorithm) {
        case BY_DEFINITION:
            resultSize = size;
            dataOUT = malloc(resultSize * sizeof(float complex));
            fftByDefinition(data, dataOUT, size);
            break;
        case COOLEY_TUKEY_RECURSIVE:
            resultSize = adjustSize(size);
            dataOUT = malloc(resultSize * sizeof(float complex));
            if(size == resultSize) {
                // Array is of power of 2
                fftCooleyTukeyRecursive(data, dataOUT, size);
            } else {
                // Array is NOT of power of 2
                float complex* adjustedData = createDataArray(size, resultSize);
                for(size_t i=0; i<size; i++) {
                    adjustedData[i] = data[i];
                }
                fftCooleyTukeyRecursive(adjustedData, dataOUT, size);
                free(adjustedData);
            }
            break;
        case COOLEY_TUKEY_ITERATIVE:
			resultSize = adjustSize(size);
            dataOUT = malloc(resultSize * sizeof(float complex));
            if(size == resultSize) {
                // Array is of power of 2
                fftCooleyTukeyIterative(data, dataOUT, size);
            } else {
                // Array is NOT of power of 2
                float complex* adjustedData = createDataArray(size, resultSize);
                for(size_t i=0; i<size; i++) {
                    adjustedData[i] = data[i];
                }
                fftCooleyTukeyIterative(adjustedData, dataOUT, size);
                free(adjustedData);
            }
            break;
        case BLUESTEIN:
            resultSize = size;
            dataOUT = malloc(resultSize * sizeof(float complex));
            size_t powerOfTwoSize = adjustSize(size * 2 + 1);
            float complex* twiddleFactors = malloc(size * sizeof(float complex));
            float complex* a = calloc(powerOfTwoSize, sizeof(float complex)); // Remember: This array must be initialized with 0
            float complex* b = calloc(powerOfTwoSize, sizeof(float complex)); // Remember: This array must be initialized with 0
            float complex* c = malloc(powerOfTwoSize * sizeof(float complex));
            float complex* d = malloc(powerOfTwoSize * sizeof(float complex));

            fftBluestein(data, dataOUT, twiddleFactors, a, b, c, d, powerOfTwoSize, size);

            free(twiddleFactors);
            free(a);
            free(b);
            free(c);
            break;
        case CHIRP_Z:
            resultSize = size;
            float complex* chirps = createDataArray(size, size);
            float complex* weights = createDataArray(size, size);
            dataOUT = malloc(resultSize * sizeof(float complex));
            fftChirpZ(data, dataOUT, chirps, weights, size);

            free(chirps);
            free(weights);
            break;
        default:
            printf("Unknown FFT algorithm! Abort...\n");
            dataOUT = NULL;
    }
    if(dataOUT == NULL) {
        return NULL;
    } else {
        struct FFTResult* result;
        result = malloc(sizeof(struct FFTResult));
        result -> data = dataOUT;
        result -> size = resultSize;
        return result;
    }
}

struct FFTResult* inverseFFT(enum FFTAlgorithm algorithm, const float complex* data, size_t size) {
	size_t resultSize;
    float complex* dataOUT;
    switch (algorithm) {
        case BY_DEFINITION:
            resultSize = size;
            dataOUT = malloc(resultSize * sizeof(float complex));
            inverseFFTByDefinition(data, dataOUT, size);
            break;
        case COOLEY_TUKEY_RECURSIVE:
            resultSize = adjustSize(size);
            dataOUT = malloc(resultSize * sizeof(float complex));
            if(size == resultSize) {
                // Array is of power of 2
                inverseFFTCooleyTukeyRecursive(data, dataOUT, size);
            } else {
                // Array is NOT of power of 2
                float complex* adjustedData = createDataArray(size, resultSize);
                for(size_t i=0; i<size; i++) {
                    adjustedData[i] = data[i];
                }
                inverseFFTCooleyTukeyRecursive(adjustedData, dataOUT, size);
                free(adjustedData);
            }
            break;
        case COOLEY_TUKEY_ITERATIVE:
            resultSize = adjustSize(size);
            dataOUT = malloc(resultSize * sizeof(float complex));
            if(size == resultSize) {
                // Array is of power of 2
                inverseFFTCooleyTukeyIterative(data, dataOUT, size);
            } else {
                // Array is NOT of power of 2
                float complex* adjustedData = createDataArray(size, resultSize);
                for(size_t i=0; i<size; i++) {
                    adjustedData[i] = data[i];
                }
                inverseFFTCooleyTukeyIterative(adjustedData, dataOUT, size);
                free(adjustedData);
            }
            break;
        case BLUESTEIN:
            resultSize = size;
            dataOUT = malloc(resultSize * sizeof(float complex));
            size_t powerOfTwoSize = adjustSize(size * 2 + 1);
            float complex* twiddleFactors = malloc(size * sizeof(float complex));
            float complex* a = calloc(powerOfTwoSize, sizeof(float complex));
            float complex* b = calloc(powerOfTwoSize, sizeof(float complex));
            float complex* c = malloc(powerOfTwoSize * sizeof(float complex));
            float complex* d = malloc(powerOfTwoSize * sizeof(float complex));

            inverseFFTBluestein(data, dataOUT, twiddleFactors, a, b, c, d, powerOfTwoSize, size);

            free(twiddleFactors);
            free(a);
            free(b);
            free(c);
            break;
        case CHIRP_Z:
            resultSize = size;
            dataOUT = malloc(size * sizeof(float complex));
            float complex* chirps = malloc(size * sizeof(float complex));
            float complex* weights = malloc(size * sizeof(float complex));

            inverseFFTChirpZ(data, dataOUT, chirps, weights, size);

            free(chirps);
            free(weights);
            break;
        default:
            printf("Unknown FFT algorithm! Abort...\n");
            dataOUT = NULL;
    }
    if(dataOUT == NULL) {
        return NULL;
    } else {
        struct FFTResult* result;
        result = malloc(sizeof(struct FFTResult));
        result -> data = dataOUT;
        result -> size = resultSize;
        return result;
    }
}