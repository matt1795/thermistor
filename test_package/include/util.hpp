// testing utilities
//
// Author: Matthew Knight
// File Name: util.hpp
// Date: 2019-09-04

#pragma once

#include "ertj0ev474j.hpp"

#include "thermistor/steinhart.hpp"

constexpr double kelvin(double celcius) {
    return celcius + Thermistor::kelvin;
}
