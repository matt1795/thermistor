// Circuit Tests
//
// Author: Matthew Knight
// File Name: circuit.hpp
// Date: 2019-09-07

#include "typical.hpp"

#include "thermistor/circuit.hpp"
#include "thermistor/ntc.hpp"

#include <gcem.hpp>
#include <gtest/gtest.h>

#include <limits>
#include <random>
#include <utility>

using TempRange = Thermistor::Range<-10, 50>;

TEST(CircuitTests, None) {
    Thermistor::Circuit::None circuit;

    std::mt19937 gen;
    std::uniform_real_distribution<double> dist{
        0.0, std::numeric_limits<double>::max()};

    for (auto i = 0; i < 10000; i++) {
        double value = dist(gen);
        EXPECT_DOUBLE_EQ(value, circuit.transform(value));
    }

    // test integration with ntc -- just make sure we  get 25C for 10000 ohms
    constexpr Thermistor::Datapoint nominal{25.0, 10000.0};
    constexpr double beta = 3950.0;
    constexpr Thermistor::Ntc<Thermistor::Range<-10, 50>, 61, double> lut{
        Thermistor::Steinhart{nominal, beta}};

    auto [expected_temp, sat] = lut.interpolate(10000);

    EXPECT_FALSE(sat);
    EXPECT_DOUBLE_EQ(expected_temp, nominal.temp);
}

TEST(CircuitTests, Adc) {
    double vref = 5.0;
    double impedance = 10000.0;

    Thermistor::Circuit::Adc<12> adc{vref};
    Thermistor::Circuit::Adc<16> adc_imp{vref, impedance};

    EXPECT_DOUBLE_EQ(adc.impedance, std::numeric_limits<double>::infinity());
    EXPECT_DOUBLE_EQ(adc_imp.impedance, impedance);

    // TODO: validate adc conversion
}

TEST(CircuitTests, HalfBridge) {
    constexpr double supply = 3.3;
    constexpr Thermistor::Circuit::HalfBridge bridge{
        Thermistor::Circuit::Adc<12>{supply}, supply, 3000.0};

    // TODO: test circuit functionality

    // TODO: test integration with ntc
    constexpr Thermistor::Ntc<Thermistor::Range<-10, 110>, 121, double,
                              std::uint16_t>
        lut{Typical::equation, bridge};
}
