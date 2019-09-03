// Steinhart lookup table
//
// Author: Matthew Knight
// File Name: steinhart.hpp
// Date: 2019-09-03

#pragma once

#include "gcem.hpp"

#include <algorithm>
#include <array>
#include <tuple>

namespace Thermistor {
    struct Datapoint {
        double temp;
        double res;
    };

    template <auto minimum, auto maximum> struct Range {
        static_assert(min < max, "min is not less than max");

        constexpr auto min = minimum;
        constexpr auto max = maximum;
    };

    // alias for readability
    template <auto minimum, auto maximum> using Beta = Range;

    double reverse_beta(double beta, double temp, Datapoint const &nominal) {
        return nominal.res *
               gcem::exp(beta / ((1.0 / temp2) - (1.0 / nominal.temp)));
    }

    struct Steinhart {
        double a, b, c;

      public:
        // calculates absolute temperature
        constexpr double calculate_temp(double res) {
            return 1 / (a + (b * gcem::log(res)) +
                        (c * gcem::pow(gcem::log(res), 3.0)));
        }

        // calculate resistance from absolute temperature
        constexpr double calculate_res(double temp) {
            double y = (1.0 / (2.0 * c)) * (a - (1.0 / temp));
            double x =
                gcem::sqrt(gcem::pow(b / (3.0 * c), 3.0) + gcem::pow(y, 2.0));

            return gcem::exp(gcem::pow(x - y, 1.0 / 3.0) -
                             gcem::pow(x + y, 1.0 / 3.0));
        }
    }

    // beta model
    template <typename TempRange, auto datapoints>
    struct Ntc {
        Output table[datapoints];

        static constexpr auto delta =
            static_cast<double>(TempRange::max - TempRange::min + 1) /
            datapoints;

        template <typename SteinhartType>
        constexpr Ntc(SteinhartType &&equation) {
            for (auto i = 0; i < datapoints; i++)
                table[i] = equation.calculate_res((i * delta) + TempRange::min +
                                                  273.15);
        }
    };

    // Regular steinhart coefficients
    template <typename TempRange, auto datapoints>
    constexpr auto make_lut(double a, double b, double c) {
        return Ntc<TempRange, datapoints>{Steinhart{a, b, c}};
    }

    // Single Beta
    template <typename TempRange, auto datapoints>
    constexpr auto make_lut(Datapoint const &nominal, double beta) {
        return Ntc<TempRange, datapoints>{Steinhart{
            (1.0 / nominal.temp) - ((1.0 / beta) * gcem::log(nominal.res)),
            1.0 / beta, 0.0}};
    }

    // Two Betas
    template <typename TempRange, auto datapoints, typename Beta1,
              typename Beta2>
    constexpr auto make_lut(Datapoint const &nominal, double beta1,
                            double beta2) {
        // make sure that beta is with respect to nominal
        // calculate three datapoints
        Datapoint one{Beta1::max, reverse_beta(beta1, Beta1::max, nominal)};
        Datapoint two{Beta2::max, reverse_beta(beta2, Beta2::max, nominal)};
        // calculate a b c from data points
        // create Steinhart
        // create Ntc
    }
} // namespace Thermistor
