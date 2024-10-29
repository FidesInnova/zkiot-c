#ifndef ZKP_H
#define ZKP_H

#include <vector>
#include <cstdint>
#include <cmath>

#include <string>


class ZKP {
public:
    void setup(int64_t g, int64_t d, int64_t l, int64_t p, std::vector<int64_t>& pp);
    void createMat(const std::string& filename, std::vector<std::vector<int64_t>>& A, std::vector<std::vector<int64_t>>& B, std::vector<std::vector<int64_t>>& C, int64_t p);
    void AHPverify(int64_t mod);
};

#endif // ZKP_H
