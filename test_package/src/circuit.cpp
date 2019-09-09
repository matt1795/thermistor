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
        std::numeric_limits<double>::min(), std::numeric_limits<double>::max()};

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

// factory for test function for ADC
template <typename AdcType>
auto create_adc_test(AdcType const& adc) {
    return [&]() {
        std::mt19937 gen;
        std::uniform_real_distribution<double> dist{-0.5, adc.vref + 0.5};
        std::uint32_t max_value = (1 << adc.resolution) - 1;
        double delta_voltage = adc.vref / max_value;

        for (auto i = 0; i < 10000; i++) {
            double value = dist(gen);
            if (value < delta_voltage)
                EXPECT_EQ(0, adc.convert(value));
            else if (value >= adc.vref)
                EXPECT_EQ(max_value, adc.convert(value));
            else
                EXPECT_EQ(static_cast<std::uint32_t>(value / delta_voltage),
                          adc.convert(value));
        }
    };
}

TEST(CircuitTests, Adc) {
    constexpr double vref = 5.0;
    constexpr double impedance = 10000.0;

    constexpr Thermistor::Circuit::Adc<12> adc{vref};
    constexpr Thermistor::Circuit::Adc<16> adc_imp{vref, impedance};

    EXPECT_DOUBLE_EQ(adc.impedance, std::numeric_limits<double>::infinity());
    EXPECT_DOUBLE_EQ(adc_imp.impedance, impedance);

    create_adc_test(adc)();
    create_adc_test(adc_imp)();
}

double bridge_calculate(double supply, double r1, double r2,
                        std::uint32_t resolution) {
    return gcem::floor(((supply * r2) / (r1 + r2)) *
                       (((1 << resolution) - 1) / supply));
}

TEST(CircuitTests, HalfBridge) {
    constexpr double supply = 3.3;
    constexpr double nominal_res = 3000.0;
    constexpr auto resolution = 12;
    constexpr auto max_value = (1 << resolution) - 1;
    constexpr Thermistor::Circuit::HalfBridge bridge{
        Thermistor::Circuit::Adc<resolution>{supply}, supply, nominal_res};

    constexpr double adc_impedance = 50000.0;
    constexpr Thermistor::Circuit::HalfBridge unbuffered_bridge{
        Thermistor::Circuit::Adc<resolution>{supply, adc_impedance}, supply,
        nominal_res};

    // TODO: test circuit functionality
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist{0.0, 10000.0};

    for (auto i = 0; i < 10000; i++) {
        double r2 = dist(gen);
        EXPECT_DOUBLE_EQ(bridge_calculate(supply, bridge.r1, r2, resolution),
                         bridge.transform(r2));

        r2 = (r2 * adc_impedance) / (r2 + adc_impedance);
        EXPECT_DOUBLE_EQ(bridge_calculate(supply, bridge.r1, r2, resolution),
                         bridge.transform(r2));
    }

    // test that it constexpr compiles with Ntc
    constexpr Thermistor::Ntc<Thermistor::Range<-10, 110>, 121, double,
                              std::uint16_t>
        lut{Typical::equation, bridge};
}
