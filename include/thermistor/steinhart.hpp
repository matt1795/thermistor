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

    // beta model
    template <typename TempRange, auto datapoints, typename Temp,
              typename Res = std::uint32_t,
              typename = std::enable_if_t<std::is_signed_v<Temp>>>
    struct Ntc {
        std::array<Res, datapoints> table{};

        static constexpr auto size = datapoints;
        static constexpr auto delta =
            static_cast<double>(TempRange::max - TempRange::min) /
            (datapoints - 1);

        template <typename IndexType>
        constexpr Temp index_to_temp(IndexType i) const {
            return static_cast<Temp>(i) * delta + TempRange::min;
        }

        template <typename Iterator>
        constexpr Temp iterator_to_temp(Iterator const& it) const {
            return index_to_temp(std::distance(it, table.rend()) - 1);
        }

        // outputs interpolated temperature and whether it is a saturated value
        std::pair<Temp, bool> interpolate(Res const& res) const {
            auto it = std::lower_bound(table.rbegin(), table.rend(), res);

            // saturate the value if out of bounds
            if (it == table.rbegin()) {
                // handle case where reading is on edge of max temp
                Temp temp = iterator_to_temp(table.rbegin());
                if (res == *table.rbegin())
                    return std::make_pair(temp, false);
                else
                    return std::make_pair(temp, true);
            } else if (it == table.rend()) {
                return std::make_pair(iterator_to_temp(std::prev(table.rend())),
                                      true);
            } else {
                // interpolate
                Temp x1 = iterator_to_temp(it);
                Temp x2 = iterator_to_temp(std::prev(it));
                Res y1 = *it;
                Res y2 = *std::prev(it);

                return std::make_pair(x1 + ((y1 - res) * (x2 - x1) / (y1 - y2)),
                                      false);
            }
        }

        constexpr Ntc(Steinhart const& equation) {
            for (auto i = 0; i < datapoints; i++) {
                double res = equation.calculate_res(
                    static_cast<double>(i * delta) + TempRange::min + kelvin);

                if constexpr (std::is_integral_v<Res>)
                    table[i] = gcem::round(res);
                else
                    table[i] = res;
            }

            if (!descending(std::begin(table), std::end(table))) {
                if (over_sampled(std::begin(table), std::end(table)))
                    throw std::logic_error(
                        "the thermistor transfer function "
                        "is over sampled and not able to distinguish between "
                        "some temperatures (decrease number of datapoints)");
                throw std::logic_error(
                    "table values must be in descending order");
            }
        }
    };

    // Regular steinhart coefficients
    template <typename TempRange, auto datapoints, typename Temp,
              typename Res = std::uint32_t>
    constexpr auto make_lut(double a, double b, double c) {
        return Ntc<TempRange, datapoints, Temp, Res>{Steinhart{a, b, c}};
    }

    // Single Beta
    template <typename TempRange, auto datapoints, typename Temp,
              typename Res = std::uint32_t>
    constexpr auto make_lut(Datapoint const& nominal, double beta) {
        double b = 1.0 / beta;
        double a =
            (1.0 / (nominal.temp + kelvin)) - (b * gcem::log(nominal.res));
        double c = 0.0;

        return make_lut<TempRange, datapoints, Temp, Res>(a, b, c);
    }

    // Two Betas
    template <typename TempRange, auto datapoints, typename Temp,
              typename Res = std::uint32_t>
    constexpr auto make_lut(Datapoint const& nominal, BetaPoint const& b1,
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
        double c = ((p2 - p1) / (l2 - l1)) * (1.0 / (l0 + l1 + l2));
        double b = p1 - (c * ((l0 * l0) + (l0 * l1) + (l1 * l1)));
        double a = y0 - l0 * (b + (c * (l0 * l0)));

        return make_lut<TempRange, datapoints, Temp, Res>(a, b, c);
    }
} // namespace Thermistor
