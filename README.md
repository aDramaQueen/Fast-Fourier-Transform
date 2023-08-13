# Fast-Fourier-Transform
Small collection of FFT algorithms that do not create additional arrays on the stack.

This small project originates from an internship in my studies:

>- The task was to program an FFT algorithm that only minimally fills the stack with new elements (at max some int/flout values, no arrays at all).
>- Everything should be preinitialized over the heap.
>- Also, the code should run on a 32 bit FPGA CPU.

This forced me to write the FFT by myself. It's not optimized in any other way, than the FFT algorithm it self, meaning there is probably room to accelerate this code by some clever C magic.

I tried to be as readable as possible. Meaning, lots of documentation & comments.

## Implemented Algorithms

- [Discrete Fourier Transformation (**DFT**)](https://en.wikipedia.org/wiki/Discrete_Fourier_transform)
- [FFT Cooley-Tukey (recursive - radix-2 decimation-in-time)](https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm)
- [FFT Cooley-Tukey (iterative - radix-2 decimation-in-time)](https://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm)
- [Bluestein](https://en.wikipedia.org/wiki/Chirp_Z-transform#Bluestein's_algorithm)
- [Chirp-Z](https://en.wikipedia.org/wiki/Chirp_Z-transform)

## Excerpt from test results

These tests ran on a machine with following stats:
- Intel(R) Core(TM) i5-5675C CPU @ 3.10GHz
- 16 GB Ram

```
================================================================
Duration (in sec.) for 4096 samples
Time for DFT:                                      1.8980000019s
Time for Cooley-Tukey (recursive):                 0.0040000002s
Time for Cooley-Tukey (iterative):                 0.0010000000s
Time for Bluestein:                                0.0140000004s
================================================================
```

# Requirements

- C-Compiler with C11 (or higher) standard
- [CMake](https://cmake.org/) v.3.5.0 (or higher)

# TODO

The Chirp-Z FFT does currently **not** work, since I didn't understand how the algorithm works...