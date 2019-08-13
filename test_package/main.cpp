#include "thermistor/ntc.hpp"
#include "thermistor/circuits.hpp"

#include <iostream>
#include <ratio>

constexpr auto ntc_lookup = Thermistor::Ntc<25,    // nominal temp Celcius
                                            10000, // nominal resistance
                                            25,    // beta temp first
                                            50,    // beta temp second
                                            3950,  // beta value
                                            -50,   // min temp Celcius
                                            150,   // max temp Celcius
                                            201    // number of datapoints
                                            >{
											Thermistor::Circuits::VoltageDividerBottom<5, std::ratio<1, 1>, 10000, 12, 5, std::ratio<1, 1>>::transform()
											};

int main() {
	std::uint32_t res_values[] = { 9783, 1000000, 150 };

	double i = -50.0;
	for (auto& res : ntc_lookup.table) {
		std::cout << i << ", " << res << std::endl;
		i += 1.0;
	}

	for (auto& resistance : res_values) {
		auto [temp, saturated] = ntc_lookup.interpolate(resistance);
		std::cout << "temp: " << temp << ", saturated: " << saturated << std::endl;
	}
}
