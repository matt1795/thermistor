// Steinhart lookup table
//
// Author: Matthew Knight
// File Name: steinhart.hpp
// Date: 2019-09-03

#pragma once

#include "util.hpp"

#include "gcem.hpp"

#include <algorithm>
#include <array>
#include <tuple>

namespace Thermistor {
    struct Steinhart {
        double a, b, c;

        // calculates absolute temperature
        constexpr double calculate_temp(double res) const {
            if (res <= 0.0)
                throw std::runtime_error(
                    "cannot have negative or zero resistance");
            return 1 / (a + (b * gcem::log(res)) +
                        (c * gcem::pow(gcem::log(res), 3.0)));
        }

        // calculate resistance from absolute temperature
        constexpr double calculate_res(double temp) const {
            if (temp <= 0.0)
                throw std::runtime_error(
                    "cannot have negative or zero absolute temperature");

            // if c was precisely set to zero -- single beta case
            if (c == 0.0) {
                return gcem::exp(((1.0 / temp) - a) / b);
            } else {
                double y = (1.0 / (2.0 * c)) * (a - (1.0 / temp));
                double x = gcem::sqrt(gcem::pow(b / (3.0 * c), 3.0) + (y * y));

                return gcem::exp(gcem::pow(x - y, 1.0 / 3.0) -
                                 gcem::pow(x + y, 1.0 / 3.0));
            }
        }
    };
} // namespace Thermistor
