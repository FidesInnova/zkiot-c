#ifndef ZKP_H
#define ZKP_H

#include <vector>
#include <cstdint>
#include <cmath>

class ZKP {
public:
    void setup(uint64_t g, uint64_t d, uint64_t l, uint64_t p, std::vector<uint64_t>& pp);
};

#endif // ZKP_H
