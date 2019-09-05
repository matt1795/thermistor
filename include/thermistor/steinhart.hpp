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

namespace {}

namespace Thermistor {
    constexpr auto kelvin = 273.15;

    struct Datapoint {
        double temp;
        double res;
    };

    struct BetaPoint {
        double temp;
        double beta;
    };

    template <auto minimum, auto maximum>
    struct Range {
        static_assert(minimum < maximum, "min is not less than max");

        static constexpr auto min = minimum;
        static constexpr auto max = maximum;
    };

    constexpr double reverse_beta(BetaPoint const& point,
                                  Datapoint const& nominal) {
        return nominal.res *
               gcem::exp(point.beta * ((1.0 / (point.temp + kelvin)) -
                                       (1.0 / (nominal.temp + kelvin))));
    }

    struct Steinhart {
        double a, b, c;

        // check if a value is inf or -inf
        constexpr bool eval_infinity(double val) const {
            double inf = std::numeric_limits<double>::infinity();
            return gcem::abs(val) == inf;
        }

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
                double m = 1.0 / (2.0 * c);
                double n = b / (3.0 * c);
                double y = m * (a - (1.0 / temp));
                double x = gcem::sqrt(gcem::pow(n, 3.0) + (y * y));

                return gcem::exp(gcem::pow(x - y, 1.0 / 3.0) -
                                 gcem::pow(x + y, 1.0 / 3.0));
            }
        }
    };

    // beta model
    template <typename TempRange, auto datapoints, typename Output = float>
    struct Ntc {
        std::array<Output, datapoints> table{};

        static constexpr auto size = datapoints;
        static constexpr auto delta =
            static_cast<double>(TempRange::max - TempRange::min) /
            (datapoints - 1);

        constexpr Ntc(Steinhart const& equation) {
            for (auto i = 0; i < datapoints; i++)
                table[i] = equation.calculate_res(
                    static_cast<double>(i * delta) + TempRange::min + kelvin);
        }
    };

    // Regular steinhart coefficients
    template <typename TempRange, auto datapoints, typename Output = float>
    constexpr auto make_lut(double a, double b, double c) {
        return Ntc<TempRange, datapoints, Output>{Steinhart{a, b, c}};
    }

    // Single Beta
    template <typename TempRange, auto datapoints, typename Output = float>
    constexpr auto make_lut(Datapoint const& nominal, double beta) {
        auto b = 1.0 / beta;
        auto a = (1.0 / (nominal.temp + kelvin)) - (b * gcem::log(nominal.res));
        auto c = 0.0;

        return make_lut<TempRange, datapoints, Output>(a, b, c);
    }

    // Two Betas
    template <typename TempRange, auto datapoints, typename Output = float>
    constexpr auto make_lut(Datapoint const& nominal, BetaPoint const& b1,
                            BetaPoint const& b2) {
        // absolute temperatures
        auto t0 = nominal.temp + kelvin;
        auto t1 = b1.temp + kelvin;
        auto t2 = b2.temp + kelvin;

        // reverse calculate resistances
        auto r1 = reverse_beta(b1, nominal);
        auto r2 = reverse_beta(b2, nominal);

        // intermediate calculations
        auto l0 = gcem::log(nominal.res);
        auto l1 = gcem::log(r1);
        auto l2 = gcem::log(r2);

        auto y0 = 1.0 / t0;
        auto y1 = 1.0 / t1;
        auto y2 = 1.0 / t2;

        auto p1 = (y1 - y0) / (l1 - l0);
        auto p2 = (y2 - y0) / (l2 - l0);

        // Steinhart coefficients
        auto c = ((p2 - p1) / (l2 - l1)) * (1.0 / (l0 + l1 + l2));
        auto b = p1 - (c * ((l0 * l0) + (l0 * l1) + (l1 * l1)));
        auto a = y0 - l0 * (b + (c * (l0 * l0)));

        return make_lut<TempRange, datapoints, Output>(a, b, c);
    }
} // namespace Thermistor
