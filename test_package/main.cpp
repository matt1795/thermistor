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
                                            251    // number of datapoints
                                            >{
											Thermistor::Circuits::VoltageDividerBottom<5, std::ratio<1, 1>, 10000, 12, 5, std::ratio<1, 1>>::transform()
											};

int main() {
	std::uint32_t adc_values[] = { 2000, 1000000, 40 };

	double i = -50.0;
	std::cout << "ADC values for each temperature:" << std::endl;
	for (auto& adc_value : ntc_lookup.table) {
		std::cout << i << ", " << adc_value << std::endl;
		i += ntc_lookup.delta;
	}

	std::cout << std::endl << "saturation output:" << std::endl;
	for (auto& adc_value : adc_values) {
		auto [temp, saturated] = ntc_lookup.interpolate(adc_value);
		std::cout << "value: " << adc_value << ", temp: " << temp << ", saturated: " << saturated << std::endl;
	}

	std::cout << std::endl << "delta: " << ntc_lookup.delta << std::endl;
}
