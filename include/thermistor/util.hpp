// Thermistor utility functions
//
// Author: Matthew Knight
// File Name: util.hpp
// Date: 2019-08-12

#pragma once

#include <iterator>

namespace Thermistor {
	// constexpr range checker. predicate is used to compare every element and
	// its predesesor
	template <typename Iterator, typename Predicate>
	constexpr bool all_of(Iterator first, Iterator last, Predicate p) {
		first++;
		for (;first != last; ++first)
			if (!p(*first, *std::prev(first)))
				return false;

		return true;
	}

	template <typename Iterator, typename Predicate>
	constexpr bool any_of(Iterator first, Iterator last, Predicate p) {
		first++;
		for (;first != last; ++first)
			if (p(*first, *std::prev(first)))
				return true;

		return false;
	}

	// checks if range is in ascending order
	template <typename Iterator>
	constexpr bool ascending(Iterator first, Iterator last) {
		return all_of(first, last, [](auto& current, auto& previous) {
			return current > previous;
		});
	}

	// checks if range is in descending order
	template <typename Iterator>
	constexpr bool descending(Iterator first, Iterator last) {
		return all_of(first, last, [](auto& current, auto& previous) {
			return current < previous;
		});
	}

	// checks to see if any values are equal
	template <typename Iterator>
	constexpr bool over_sampled(Iterator first, Iterator last) {
		return any_of(first, last, [](auto& current, auto& previous) {
			return current == previous;
		});
	}
}
