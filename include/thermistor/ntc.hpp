// NTC Thermistor Class
//
// Author: Matthew Knight
// File Name: ntc.hpp
// Date: 2019-08-08

#pragma once

#include "gcem.hpp"

#include <algorithm>
#include <array>
#include <tuple>

namespace Thermistor {
    // beta model
    template <auto nominal_temp, auto nominal_res, auto beta_first,
              auto beta_second, auto beta, auto min_temp, auto max_temp,
              auto datapoints, typename Output = float,
              typename Input = std::uint32_t>
    struct Ntc {
        Output table[datapoints];

        template <typename T>
        static constexpr double c_to_k(T const& celcius) {
            return 273.15 + celcius;
        }

        template <typename T>
        static constexpr double k_to_c(T const& kelvin) {
            return static_cast<double>(kelvin) - 273.15;
        }

        // calculations for this thermistor
        static constexpr Output calculate(Input const& res) {
            return k_to_c(1.0 / ((1.0 / c_to_k(nominal_temp)) -
                                 (gcem::log(nominal_res / res) / beta)));
        }

        static constexpr Input reverse_calc(Output const& temp) {
            return static_cast<double>(nominal_res) *
                   gcem::exp(
                       -1.0 * beta *
                       ((1.0 / c_to_k(nominal_temp)) - (1.0 / c_to_k(temp))));
        }

        static auto constexpr var = -1.0 * beta;
        static auto constexpr max_temp_res = reverse_calc(max_temp);
        static auto constexpr min_temp_res = reverse_calc(min_temp);
        static auto constexpr delta =
            gcem::abs(static_cast<double>(max_temp_res) - min_temp_res) /
            datapoints;

        constexpr Ntc() {
            if constexpr (nominal_temp == beta_first) {
                for (auto i = 0; i < datapoints; i++)
                    table[i] = calculate((i * delta) + max_temp_res);
            }
        }

        // returns temperature and whether the value is saturated
        constexpr Output lookup(Input const& input) const { return 0.0; }

        constexpr auto interpolate() const { return 0.0; }
    };
} // namespace Thermistor
