// Thermistor Tests
//
// Author: Matthew Knight
// File Name: thermistor.cpp
// Date: 2019-09-04

// Tests for this library

#include <gcem.hpp>
#include <gtest/gtest.h>

#include "ertj0ev474j.hpp"
#include "typical.hpp"

#include "thermistor/steinhart.hpp"

#include <array>
#include <tuple>

TEST(ThermistorTests, SteinhartTest) {
    for (auto& datapoint : Typical::data) {
        // use foat comparison to handle tiny rounding error
        EXPECT_FLOAT_EQ(Typical::equation.calculate_res(datapoint.first),
                        datapoint.second);
        EXPECT_FLOAT_EQ(Typical::equation.calculate_temp(datapoint.second),
                        datapoint.first);
    }
}

TEST(ThermistorTests, SteinhartLookupTest) {
    using TempRange = Thermistor::Range<-10, 50>;
    constexpr auto lut =
        Thermistor::Ntc<TempRange, Typical::data.size()>(Typical::equation);

    for (auto i = 0; i < lut.table.size(); i++)
        EXPECT_FLOAT_EQ(lut.table[i], Typical::data[i].second);
}

TEST(ThermistorTests, MakeLutTest) {
    using TempRange = Thermistor::Range<-10, 50>;
    constexpr auto lut = Thermistor::make_lut<TempRange, Typical::data.size()>(
        Typical::equation.a, Typical::equation.b, Typical::equation.c);

    for (auto i = 0; i < lut.table.size(); i++) {
        EXPECT_FLOAT_EQ(lut.table[i], Typical::data[i].second);
    }
}

TEST(ThermistorTests, SingleBetaTest) {
    // create lut for single beta, double check table values against equation
    using TempRange = Thermistor::Range<-55, 105>;
    constexpr auto beta = 3892.0;
    constexpr Thermistor::Datapoint nominal{25.0, 10000.0};
    constexpr auto datapoints = 161;
    constexpr auto lut =
        Thermistor::make_lut<TempRange, datapoints>(nominal, beta);

    for (auto i = 0; i < lut.table.size(); i++) {
        auto temp = static_cast<double>(i * lut.delta) + TempRange::min +
                    Thermistor::kelvin;
        double expected_res =
            nominal.res * gcem::exp(beta * ((1.0 / temp) - 1.0 / kelvin(nominal.temp)));

        EXPECT_FLOAT_EQ(expected_res, lut.table[i]);
    }
}

template <typename LutType, typename Data>
auto mean_squared_error(LutType const& lut, Data const& table) {
    if (LutType::size != table.size())
        throw std::runtime_error("sizes do not match");

    double acc = 0.0;
    for (auto i = 0; i < table.size(); i++) {
        auto err = gcem::abs(table[i].second - lut.table[i]);
        acc += err * err;
    }

    return acc / table.size();
}

TEST(ThermistorTests, DoubleBetaTest) {
    // create lut for double beta, compare against actual table for thermistor,
    // ensure that squared error is smaller than using single betas
    using TempRange = Thermistor::Range<-40, 125>;
    constexpr auto datapoints = TempRange::max - TempRange::min + 1;
    constexpr Thermistor::Datapoint nominal{25.0, 470000.0};

    constexpr auto beta1 = 4700.0;
    constexpr auto beta2 = 4750.0;
    constexpr Thermistor::BetaPoint b1{50.0, beta1};
    constexpr Thermistor::BetaPoint b2{85.0, beta2};

    constexpr auto single1 =
        Thermistor::make_lut<TempRange, datapoints>(nominal, beta1);
    constexpr auto single2 =
        Thermistor::make_lut<TempRange, datapoints>(nominal, beta2);
    constexpr auto full =
        Thermistor::make_lut<TempRange, datapoints>(nominal, b1, b2);

    auto err_single1 = mean_squared_error(single1, ertj0ev474j::data);
    auto err_single2 = mean_squared_error(single2, ertj0ev474j::data);
    auto err_full = mean_squared_error(full, ertj0ev474j::data);

    // Test that double beta mse is less than using either single beta
    EXPECT_TRUE(err_full < err_single1 && err_full < err_single2);
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
