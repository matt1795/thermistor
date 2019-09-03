#include "thermistor/ntc.hpp"

#include <iostream>

constexpr auto ntc_lookup = Thermistor::Ntc<25,    // nominal temp Celcius
                                            10000, // nominal resistance
                                            25,    // beta temp first
                                            50,    // beta temp second
                                            3950,  // beta value
                                            -50,   // min temp Celcius
                                            150,   // max temp Celcius
                                            201    // number of datapoints
                                            >{};

constexpr auto example = Thermistor::make_lut<
    
int main() {
    std::cout << "delta: " << ntc_lookup.delta << std::endl;
    std::cout << ntc_lookup.lookup(10000) << std::endl;

    int i = 0;
    for (auto& datapoint : ntc_lookup.table)
        std::cout << ((i++ * ntc_lookup.delta) + -50) << ": " <<datapoint << std::endl;
}
