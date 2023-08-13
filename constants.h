/**
 * Application wide constant values
 * 
 * @version 1.0
 * @date 2023-07-10
 * @file fft.h
 * @author Richard Saeuberlich (richard.saeuberlich@stud.tu-darmstadt.de)
 */
#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#include <stdbool.h>

/**
 * Runs this application in Debug mode.
 * This activates some sanity tests, that may be omitted during production.
 */
const bool DEBUG = false;

/**
 * Downcast Pi to 32 bit floating point, since this application is supposed
 * to work with 32 bit floating point numbers.
 */
const float PI = (float)M_PI;