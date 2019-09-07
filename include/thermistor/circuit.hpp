// Value transformations on thermistor table to reflect circuitry
//
// Author: Matthew Knight
// File Name: circuits.hpp
// Date: 2019-08-12

#pragma once

#include "gcem.hpp"

namespace Thermistor::Circuit {
    struct None {
        static constexpr double transform(double res) {
            return res;
        }
    };
    template <auto resolution, auto vref, typename VrefRatio>
    struct Adc {
        // TODO: impedance?
        template <typename T>
        static constexpr std::uint32_t convert(T val) {
            // TODO: static assert that the resolution isn't too big
            double ratio = 1.0 / vref / VrefRatio::num * VrefRatio::den * val;
            if (ratio > 1.0)
                ratio = 1.0;
            if (ratio < 0.0)
                ratio = 0.0;

            return ratio * ((1 << resolution) - 1);
        }
    };

    // A half-bridge assumes that the thermistor is connected to ground.
    template <auto supply, typename SupplyRatio, auto top_res, auto resolution,
              auto vref, typename VrefRatio>
    struct HalfBridge {
        static constexpr auto transform() {
            return [](auto val) {
                double voltage = (static_cast<double>(supply) *
                                  SupplyRatio::num * SupplyRatio::den * val) /
                                 (top_res + val);
                return Adc<resolution, vref, VrefRatio>::convert(voltage);
            };
        }
    };
} // namespace Thermistor::Circuits
