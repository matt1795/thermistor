#include "thermistor/ntc.hpp"

#include <iostream>

constexpr auto ntc_lookup = Thermistor::Ntc<
	25,		// nominal temp Celcius
	10000,  // nominal resistance
	25, 	// beta temp first
	50,		// beta temp second
	3950,	// beta value
	-50,	// min temp Celcius
	150, 	// max temp Celcius
	201 	// number of datapoints
>{};

int main() {
    std::cout << "var: " << ntc_lookup.reverse_calc(24) << std::endl;
    std::cout << "delta: " << ntc_lookup.delta << std::endl;
    std::cout << "max_temp_res: " << ntc_lookup.max_temp_res << std::endl;
    std::cout << "min_temp_res: " << ntc_lookup.min_temp_res << std::endl;
	std::cout << ntc_lookup.lookup(10000) << std::endl;

    int i = 0;
    for (auto& datapoint : ntc_lookup.table) {
        std::cout << ((i++ * ntc_lookup.delta) + ntc_lookup.max_temp_res ) << ": " << datapoint << std::endl;


    }
}
