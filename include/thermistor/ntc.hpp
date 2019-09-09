// NTC Thermistor Class
//
// Author: Matthew Knight
// File Name: ntc.hpp
// Date: 2019-08-08

#pragma once

#include "circuit.hpp"
#include "steinhart.hpp"
#include "util.hpp"

#include "gcem.hpp"

#include <algorithm>
#include <array>
#include <tuple>

namespace Thermistor {
    template <auto minimum, auto maximum>
    struct Range {
        static_assert(minimum < maximum, "min is not less than max");

        static constexpr auto min = minimum;
        static constexpr auto max = maximum;
    };

    template <typename TempRange, auto datapoints, typename Temp,
              typename TableValue = std::uint32_t,
              typename = std::enable_if_t<std::is_signed_v<Temp>>>
    class Ntc {
        using Table = std::array<TableValue, datapoints>;
        Table table{};

      public:
        static constexpr auto delta =
            static_cast<double>(TempRange::max - TempRange::min) /
            (datapoints - 1);

        template <typename Circuit>
        constexpr Ntc(Steinhart const& equation,
                      Circuit const& circuit = Thermistor::Circuit::None{}) {
            for (auto i = 0; i < datapoints; i++) {
                double res = equation.calculate_res(
                    static_cast<double>(i * delta) + TempRange::min + kelvin);

                auto value = circuit.transform(res);

                if constexpr (std::is_integral_v<TableValue>)
                    table[i] = gcem::round(value);
                else
                    table[i] = value;
            }

            if (!descending(std::begin(table), std::end(table))) {
                if (over_sampled(std::begin(table), std::end(table)))
                    throw std::logic_error(
                        "the thermistor transfer function is over sampled "
                        "and not able to distinguish between some "
                        "temperatures (decrease number of datapoints)");
                throw std::logic_error(
                    "table values must be in descending order");
            }
        }

        constexpr Ntc(Steinhart const& equation)
            : Ntc(equation, Circuit::None{}) {}

        template <typename IndexType>
        constexpr Temp index_to_temp(IndexType i) const {
            return static_cast<Temp>(i) * delta + TempRange::min;
        }

        template <typename Iterator>
        constexpr Temp iterator_to_temp(Iterator const& it) const {
            return index_to_temp(std::distance(it, table.rend()) - 1);
        }

        constexpr auto cbegin() const noexcept { return table.cbegin(); }

        constexpr auto cend() const noexcept { return table.cend(); }

        constexpr auto size() const noexcept { return table.size(); }

        constexpr auto operator[](typename Table::size_type pos) const {
            return table[pos];
        }

        // outputs interpolated temperature and whether it is a saturated
        // value
        std::pair<Temp, bool> interpolate(TableValue const& res) const {
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
                TableValue y1 = *it;
                TableValue y2 = *std::prev(it);

                return std::make_pair(x1 + ((y1 - res) * (x2 - x1) / (y1 - y2)),
                                      false);
            }
        }
    };
} // namespace Thermistor
