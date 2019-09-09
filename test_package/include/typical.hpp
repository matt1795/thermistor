// Typical 3k thermistor table
//
// Author: Matthew Knight
// File Name: typical.hpp
// Date: 2019-09-05

#pragma once

#include "util.hpp"

#include "thermistor/ntc.hpp"

#include <array>
#include <tuple>

namespace Typical {
    constexpr Thermistor::Steinhart equation{1.4e-3, 2.37e-4, 9.9e-8};

    constexpr std::array data{std::make_pair(kelvin(-10.0), 17000.007273577237),
                              std::make_pair(kelvin(0.0), 10030.217595116217),
                              std::make_pair(kelvin(10.0), 6110.3822188932181),
                              std::make_pair(kelvin(20.0), 3833.3275927458376),
                              std::make_pair(kelvin(30.0), 2470.6001126460119),
                              std::make_pair(kelvin(40.0), 1632.3750999273645),
                              std::make_pair(kelvin(50.0), 1103.5474105037119)};
} // namespace Typical
