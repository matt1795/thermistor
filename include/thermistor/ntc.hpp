// NTC Thermistor Class
//
// Author: Matthew Knight
// File Name: ntc.hpp
// Date: 2019-08-08

#pragma once

#include "util.hpp"

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

        // calculate resistance for given temperature
        static constexpr Input calculate(Output const& temp) {
            return static_cast<double>(nominal_res) *
                   gcem::exp(
                       -1.0 * beta *
                       ((1.0 / c_to_k(nominal_temp)) - (1.0 / c_to_k(temp))));
        }

        template <typename IndexType>
        static constexpr double index_to_temp(IndexType i) {
            return static_cast<double>(i) * delta + min_temp;
        }

        template <typename Iterator>
        Output iterator_to_temp(Iterator const& it) const {
            return index_to_temp(std::distance(it, std::rend(table) - 1));
        }

        static auto constexpr delta =
            static_cast<double>(max_temp - min_temp) / (datapoints - 1);

        template <typename Func>
        constexpr Ntc(Func f) {
            static_assert(datapoints > 1, "datapoints must be greater than 1");

            if constexpr (nominal_temp == beta_first) {
                for (auto i = 0; i < datapoints; i++)
                    table[i] = f(calculate(index_to_temp(i)));
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

        constexpr Ntc()
            : Ntc([](auto val) { return val; }) {}

        // outputs interpolated temperature and whether it is a saturated value
        constexpr std::pair<Output, bool>
        interpolate(Input const& input) const {
            auto it =
                std::lower_bound(std::rbegin(table), std::rend(table), input);

            // saturate the value if out of bounds
            if (it == std::rbegin(table)) {
                return std::make_pair(iterator_to_temp(std::rbegin(table)),
                                      true);
            } else if (it == std::rend(table)) {
                return std::make_pair(
                    iterator_to_temp(std::prev(std::rend(table))), true);
            } else {
                // interpolate
                auto x1 = iterator_to_temp(it);
                auto y1 = *it;
                auto x2 = iterator_to_temp(std::prev(it));
                auto y2 = *std::prev(it);

                return std::make_pair(
                    (((input - y1) * (x2 - x1)) / (y2 - y1)) + x1, false);
            }
        }
    };
} // namespace Thermistor
