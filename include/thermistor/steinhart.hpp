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
    constexpr auto kelvin = 273.15;

    struct Datapoint {
        double temp;
        double res;
    };

    struct BetaPoint {
        double temp;
        double beta;
    };

    constexpr double reverse_beta(BetaPoint const& point,
                                  Datapoint const& nominal) {
        return nominal.res *
               gcem::exp(point.beta * ((1.0 / (point.temp + kelvin)) -
                                       (1.0 / (nominal.temp + kelvin))));
    }

    class Steinhart {
        double a{};
        double b{};
        double c{};

      public:
        constexpr Steinhart(double a, double b, double c)
            : a(a)
            , b(b)
            , c(c) {}

        constexpr Steinhart(Datapoint const& nominal, double beta)
            : a((1.0 / (nominal.temp + kelvin)) -
                ((1.0 / beta) * gcem::log(nominal.res)))
            , b(1.0 / beta)
            , c(0.0)

        {}

        constexpr Steinhart(Datapoint const& nominal, BetaPoint const& b1,
                            BetaPoint const& b2) {
            // absolute temperatures
            double t0 = nominal.temp + kelvin;
            double t1 = b1.temp + kelvin;
            double t2 = b2.temp + kelvin;

            // reverse calculate resistances
            double r1 = reverse_beta(b1, nominal);
            double r2 = reverse_beta(b2, nominal);

            // intermediate calculations
            double l0 = gcem::log(nominal.res);
            double l1 = gcem::log(r1);
            double l2 = gcem::log(r2);

            double y0 = 1.0 / t0;
            double y1 = 1.0 / t1;
            double y2 = 1.0 / t2;

            double p1 = (y1 - y0) / (l1 - l0);
            double p2 = (y2 - y0) / (l2 - l0);

            // Steinhart coefficients
            c = ((p2 - p1) / (l2 - l1)) * (1.0 / (l0 + l1 + l2));
            b = p1 - (c * ((l0 * l0) + (l0 * l1) + (l1 * l1)));
            a = y0 - l0 * (b + (c * (l0 * l0)));
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
                double y = (1.0 / (2.0 * c)) * (a - (1.0 / temp));
                double x = gcem::sqrt(gcem::pow(b / (3.0 * c), 3.0) + (y * y));

                return gcem::exp(gcem::pow(x - y, 1.0 / 3.0) -
                                 gcem::pow(x + y, 1.0 / 3.0));
            }
        }
    };
} // namespace Thermistor
