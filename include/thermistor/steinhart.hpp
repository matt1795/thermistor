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

namespace {

}

namespace Thermistor {
    constexpr auto kelvin = 273.15;

    struct Datapoint {
        const double temp;
        const double res;

        constexpr auto get_temp() const {
            return temp;
        }
    };

    struct BetaPoint {
        const double temp;
        const double beta;
    };

    template <auto minimum, auto maximum> struct Range {
        static_assert(minimum < maximum, "min is not less than max");

        static constexpr auto min = minimum;
        static constexpr auto max = maximum;
    };

    constexpr double reverse_beta(BetaPoint const &point,
                                  Datapoint const &nominal) {
        return nominal.res * gcem::exp(point.beta / ((1.0 / point.temp) -
                                                     (1.0 / nominal.temp)));
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
            double x = gcem::sqrt(gcem::pow(b / (3.0 * c), 3.0) + (y * y));

            return gcem::exp(gcem::pow(x - y, 1.0 / 3.0) -
                             gcem::pow(x + y, 1.0 / 3.0));
        }
    };

    // beta model
    template <typename TempRange, auto datapoints, typename Output = float>
    struct Ntc {
        Output table[datapoints];

        static constexpr auto delta =
            static_cast<double>(TempRange::max - TempRange::min + 1) /
            datapoints;

        template <typename SteinhartType>
        constexpr Ntc(SteinhartType &&equation) {
            for (auto i = 0; i < datapoints; i++)
                table[i] = equation.calculate_res((i * delta) + TempRange::min +
                                                  kelvin);
        }
    };

    // Regular steinhart coefficients
    template <typename TempRange, auto datapoints, typename Output>
    constexpr auto make_lut(double a, double b, double c) {
        return Ntc<TempRange, datapoints, Output>{Steinhart{a, b, c}};
    }

    // Single Beta
    template <typename TempRange, auto datapoints, typename Output>
    constexpr auto make_lut(Datapoint const &nominal, double beta) {
        return Ntc<TempRange, datapoints, Output>{Steinhart{
            (1.0 / nominal.temp) - ((1.0 / beta) * gcem::log(nominal.res)),
            1.0 / beta, 0.0}};
    }

    // Two Betas
    /*
    template <typename TempRange, auto datapoints, typename Output = float>
    constexpr auto make_lut(Datapoint const& nominal, BetaPoint const& b1,
                            BetaPoint const& b2) {

        // absolute temperatures
        constexpr auto t0 = nominal.temp + kelvin;
        constexpr auto t1 = b1.temp + kelvin;
        constexpr auto t2 = b2.temp + kelvin;

        // reverse calculate resistances
        constexpr auto r1 = reverse_beta(b1, nominal);
        constexpr auto r2 = reverse_beta(b2, nominal);

        // intermediate calculations
        constexpr auto l0 = gcem::log(nominal.res);
        constexpr auto l1 = gcem::log(r1);
        constexpr auto l2 = gcem::log(r2);

        constexpr auto y0 = 1.0 / t0;
        constexpr auto y1 = 1.0 / t1;
        constexpr auto y2 = 1.0 / t2;

        constexpr auto p1 = (y1 - y0) / (l1 - l0);
        constexpr auto p2 = (y2 - y0) / (l2 - l0);

        // Steinhart coefficients
        constexpr auto c = ((p2 - p1) / (l2 - l1)) * (1.0 / (l0 + l1 + l2));
        constexpr auto b = p1 - (c * ((l0 * l0) + (l0 * l1) + (l1 * l1)));
        constexpr auto a = y0 - l0 * (b + (c * (l0 * l0)));

        return make_lut<TempRange, datapoints, Output>(a, b, c);
    }
    */
} // namespace Thermistor
