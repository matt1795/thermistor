// testing utilities
//
// Author: Matthew Knight
// File Name: util.hpp
// Date: 2019-09-04

#pragma once

#include "thermistor/steinhart.hpp"

constexpr double kelvin(double celcius) {
    return celcius + Thermistor::kelvin;
}
