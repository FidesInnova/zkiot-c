#include "zkp.h"

uint64_t power(uint64_t base, uint64_t exponent, uint64_t modulus) {
    uint64_t result = 1;
    base = base % modulus;
    while (exponent > 0) {
        if (exponent % 2 == 1) {
            result = (result * base) % modulus;
        }
        exponent = exponent >> 1;
        base = (base * base) % modulus;
    }
    return result;
}

void ZKP::setup(uint64_t g, uint64_t d, uint64_t l, uint64_t p, std::vector<uint64_t>& pp) {
    pp.clear();
    uint64_t pMinusOne = p - 1;
    uint64_t current_exponent = 1;
    uint64_t newPower = d % pMinusOne;
    uint64_t res = 0;

    res = power(g, newPower, p);
    pp.push_back(res);
    for (uint64_t i = 1; i < l; ++i) {
        res = power(res, newPower, p);
        pp.push_back(res);
    }
    //pp=(2,66,83,91,96,24,2,66,83)
}
