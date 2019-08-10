// NTC Thermistor Class
//
// Author: Matthew Knight
// File Name: ntc.hpp
// Date: 2019-08-08

#pragma once

#include "gcem.hpp"

#include <array>
#include <algorithm>
#include <tuple>

namespace Thermistor {
	// beta model
	template <auto nominal_temp, auto nominal_res, auto beta_first, auto beta_second, auto beta, 
		auto min_temp, auto max_temp, auto datapoints,
		typename Output = float, typename Input = std::uint32_t>
	struct Ntc {
		Output table[datapoints];
	
		// calculations for this thermistor
		static constexpr Output calculate(Input const& res) {
			return 1.0 / ((1.0 / nominal_temp) - gcem::log(static_cast<double>(nominal_res) / res));
		}
		
		static constexpr Input reverse_calc(Output const& temp) {
			return nominal_res * gcem::exp(-1.0 * beta * ((1.0  / nominal_temp) - (1.0 / temp)));
		}
	
		constexpr Ntc() {
			/*
			if constexpr (nominal_temp == beta_first) {
				auto max_temp_res = reverse_calc(max_temp);
				auto min_temp_res = reverse_calc(min_temp);
				auto delta = gcem::abs(static_cast<double>(max_temp_res) - min_temp_res) / datapoints;
				
				for (auto i = 0; i < datapoints; i++)
					table[i] = calculate((i * delta) + min_temp_res);
			}
			*/
		}
		
		// returns temperature and whether the value is saturated
		constexpr Output lookup(Input const& input) const {
			return 0.0;
		}

		constexpr auto interpolate() {
		
		}
	};
}
