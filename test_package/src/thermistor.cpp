// Thermistor Tests
//
// Author: Matthew Knight
// File Name: thermistor.cpp
// Date: 2019-09-04

// Tests for this library

#include <gcem.hpp>
#include <gtest/gtest.h>

#include "thermistor/steinhart.hpp"

#include <array>
#include <tuple>

constexpr double kelvin(double c) { return c + Thermistor::kelvin; }

TEST(ThermistorTests, SteinhartTest) {
    // table of values to test against
    std::array data{std::make_pair(kelvin(-10.0), 17000.007273577237),
                    std::make_pair(kelvin(0.0), 10030.217595116217),
                    std::make_pair(kelvin(10.0), 6110.3822188932181),
                    std::make_pair(kelvin(20.0), 3833.3275927458376),
                    std::make_pair(kelvin(30.0), 2470.6001126460119),
                    std::make_pair(kelvin(40.0), 1632.3750999273645),
                    std::make_pair(kelvin(50.0), 1103.5474105037119)};

    // coefficients from typical 3k @ 25C
    Thermistor::Steinhart equation{1.4e-3, 2.37e-4, 9.9e-8};

    for (auto &datapoint : data) {
        // use foat comparison to handle tiny rounding error
        EXPECT_FLOAT_EQ(equation.calculate_res(datapoint.first),
                        datapoint.second);
        EXPECT_FLOAT_EQ(equation.calculate_temp(datapoint.second),
                        datapoint.first);
    }
}

TEST(ThermistorTests, SingeBetaTest) {
    // create lut for single beta, double check table values against equation
    EXPECT_TRUE(false);
}

TEST(ThermistorTests, DoubleBetaTest) {
    // create lut for double beta, compare against actual table for thermistor,
    // ensure that error is small
    EXPECT_TRUE(false);
}

TEST(ThermistorTests, InterpolationTest) {
    // make sure that the correct points are selected for interpolation
    // make sure that the calculation is correct
    EXPECT_TRUE(false);
}

TEST(ThermistorTests, SaturationTest) {
    // Ensure that input values outside of
    EXPECT_TRUE(false);
}
