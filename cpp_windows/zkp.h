#ifndef ZKP_H
#define ZKP_H

#include <vector>
#include <cstdint>
#include <cmath>

#include <string>


class ZKP {
public:
    void setup(uint64_t g, uint64_t d, uint64_t l, uint64_t p, std::vector<uint64_t>& pp);
    void createMat(const std::string& filename, std::vector<std::vector<int>>& A, std::vector<std::vector<int>>& B, std::vector<std::vector<int>>& C, uint64_t p);
};

#endif // ZKP_H
