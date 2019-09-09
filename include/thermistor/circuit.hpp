// Value transformations on thermistor table to reflect circuitry
//
// Author: Matthew Knight
// File Name: circuits.hpp
// Date: 2019-08-12

#pragma once

#include "gcem.hpp"

namespace Thermistor::Circuit {
    struct None {
        constexpr double transform(double res) const { return res; }
    };

    template <auto resolution>
    struct Adc {
        double vref;
        double impedance;

        constexpr Adc(double vref, double impedance =
                                       std::numeric_limits<double>::infinity())
            : vref(vref)
            , impedance(impedance) {

            if (vref <= 0.0)
                throw std::runtime_error("vref must be greater than zero");

            if (impedance <= 0.0)
                throw std::runtime_error("impedance must be greater than zero");
        }

        constexpr double convert(double voltage) const {
            double ratio = voltage / vref;
            if (ratio > 1.0)
                ratio = 1.0;
            else if (ratio < 0.0)
                ratio = 0.0;

            return gcem::floor(ratio * ((1 << resolution) - 1));
        }
    };

    // A half-bridge assumes that the thermistor is connected to ground.
    template <typename AdcType>
    struct HalfBridge {
        AdcType adc;
        double supply;
        double r1;

        constexpr HalfBridge(AdcType const& adc, double supply, double res)
            : adc(adc)
            , supply(supply)
            , r1(res) {}

        constexpr double transform(double res) const {
            double r2 =
                (adc.impedance == std::numeric_limits<double>::infinity())
                    ? res
                    : ((adc.impedance * res) / (adc.impedance + res));

            return adc.convert((supply * r2) / (r1 + r2));
        }
    };
} // namespace Thermistor::Circuit
