// Thermistor Tests
//
// Author: Matthew Knight
// File Name: thermistor.cpp // Date: 2019-09-04

// Tests for this library

#include <gcem.hpp>
#include <gtest/gtest.h>

#include "ertj0ev474j.hpp"
#include "typical.hpp"

#include "thermistor/ntc.hpp"

#include <array>
#include <limits>
#include <random>
#include <tuple>

template <typename LutType, typename Data>
auto mean_squared_error(LutType const& lut, Data const& table) {
    if (lut.size() != table.size())
        throw std::runtime_error("sizes do not match");

    double acc = 0.0;
    for (auto i = 0; i < table.size(); i++) {
        auto err = gcem::abs(table[i].second - lut[i]);
        acc += err * err;
    }

    return acc / table.size();
}

TEST(NtcTests, SteinhartTest) {
    for (auto& datapoint : Typical::data) {
        // use foat comparison to handle tiny rounding error
        EXPECT_FLOAT_EQ(Typical::equation.calculate_res(datapoint.first),
                        datapoint.second);
        EXPECT_FLOAT_EQ(Typical::equation.calculate_temp(datapoint.second),
                        datapoint.first);
    }
}

TEST(NtcTests, SteinhartLookupTest) {
    using TempRange = Thermistor::Range<-10, 50>;
    constexpr Thermistor::Ntc<TempRange, Typical::data.size(), double, double>
        lut{Typical::equation};

    constexpr Thermistor::Ntc<TempRange, Typical::data.size(), double,
                              std::uint32_t>
        lut_integral{Typical::equation};

    for (auto i = 0; i < lut.size(); i++)
        EXPECT_DOUBLE_EQ(lut[i], Typical::data[i].second);

    for (auto i = 0; i < lut_integral.size(); i++)
        EXPECT_EQ(lut_integral[i], gcem::round(Typical::data[i].second));
}

TEST(NtcTests, MakeLutTest) {
    using TempRange = Thermistor::Range<-10, 50>;
    constexpr Thermistor::Ntc<TempRange, Typical::data.size(), double, double>
        lut{Typical::equation};

    for (auto i = 0; i < lut.size(); i++) {
        EXPECT_DOUBLE_EQ(lut[i], Typical::data[i].second);
    }
}

TEST(NtcTests, SingleBetaTest) {
    // create lut for single beta, double check table values against equation
    using TempRange = Thermistor::Range<-55, 105>;
    constexpr double beta = 3892.0;
    constexpr Thermistor::Datapoint nominal{25.0, 10000.0};
    constexpr std::uint32_t datapoints = 161;
    constexpr Thermistor::Ntc<TempRange, datapoints, double, double> lut{
        Thermistor::Steinhart{nominal, beta}};

    for (auto i = 0; i < lut.size(); i++) {
        auto temp = static_cast<double>(i * lut.delta) + TempRange::min +
                    Thermistor::kelvin;
        double expected_res =
            nominal.res *
            gcem::exp(beta * ((1.0 / temp) - 1.0 / kelvin(nominal.temp)));

        EXPECT_FLOAT_EQ(expected_res, lut[i]);
    }
}

TEST(NtcTests, DoubleBetaTest) {
    // create lut for double beta, compare against actual table for thermistor,
    // ensure that squared error is smaller than using single betas
    using TempRange = Thermistor::Range<-40, 125>;
    constexpr auto datapoints = TempRange::max - TempRange::min + 1;
    constexpr Thermistor::Datapoint nominal{25.0, 470000.0};

    constexpr double beta1 = 4700.0;
    constexpr double beta2 = 4750.0;
    constexpr Thermistor::BetaPoint b1{50.0, beta1};
    constexpr Thermistor::BetaPoint b2{85.0, beta2};

    constexpr Thermistor::Ntc<TempRange, datapoints, float, float> single1{
        Thermistor::Steinhart{nominal, beta1}};

    constexpr Thermistor::Ntc<TempRange, datapoints, float, float> single2{
        Thermistor::Steinhart{nominal, beta2}};

    constexpr Thermistor::Ntc<TempRange, datapoints, float, float> full{
        Thermistor::Steinhart{nominal, b1, b2}};

    double err_single1 = mean_squared_error(single1, ertj0ev474j::data);
    double err_single2 = mean_squared_error(single2, ertj0ev474j::data);
    double err_full = mean_squared_error(full, ertj0ev474j::data);

    // Test that double beta mse is less than using either single beta
    EXPECT_TRUE(err_full < err_single1 && err_full < err_single2);
}

TEST(NtcTests, InterpolationTest) {
    using TempRange = Thermistor::Range<-10, 10>;
    constexpr Thermistor::Ntc<TempRange, 21, double> lut{Typical::equation};

    std::mt19937 gen;
    for (auto it = lut.begin(); it != std::prev(lut.end()); ++it) {
        std::uniform_int_distribution dist(*std::next(it), *it);

        double min_temp = std::distance(lut.begin(), it) + TempRange::min;
        double max_temp = min_temp + 1.0;

        for (auto i = 0; i < 100; i++) {
            std::uint32_t res = dist(gen);
            auto [temp, sat] = lut.interpolate(res);
            EXPECT_FALSE(sat);

            double ratio =
                static_cast<double>(*it - res) / (*it - *std::next(it));
            double expected_temp = min_temp + (ratio * (max_temp - min_temp));
            EXPECT_DOUBLE_EQ(expected_temp, temp);
        }
    }
}

// Tests the end points and just past them to catch any boundary cases
TEST(NtcTests, SaturationTest) {
    using TempRange = Thermistor::Range<0, 10>;
    constexpr Thermistor::Ntc<TempRange, 11, double> lut{Typical::equation};

    std::uint32_t max = lut[0];
    std::uint32_t min = lut[10];
    std::uint32_t max_outside = max + 10;
    std::uint32_t min_outside = min - 10;

    auto [max_temp, max_sat] = lut.interpolate(max);
    auto [min_temp, min_sat] = lut.interpolate(min);
    auto [max_temp_outside, max_sat_outside] = lut.interpolate(max_outside);
    auto [min_temp_outside, min_sat_outside] = lut.interpolate(min_outside);

    EXPECT_FALSE(max_sat);
    EXPECT_FALSE(min_sat);
    EXPECT_TRUE(max_sat_outside);
    EXPECT_TRUE(min_sat_outside);

    EXPECT_DOUBLE_EQ(max_temp, static_cast<double>(TempRange::min));
    EXPECT_DOUBLE_EQ(min_temp, static_cast<double>(TempRange::max));
    EXPECT_DOUBLE_EQ(max_temp_outside, static_cast<double>(TempRange::min));
    EXPECT_DOUBLE_EQ(min_temp_outside, static_cast<double>(TempRange::max));
}
