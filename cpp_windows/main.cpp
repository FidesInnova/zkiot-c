#include <iostream>
#include <vector>
#include "zkp.h"

int main() {
    uint64_t p = 181;
    uint64_t g = 2;
    uint64_t d = 111213119;
    uint64_t l = 8;
    uint64_t n = 5;
    std::vector<uint64_t> pp;

    ZKP lib;
    lib.setup(g, d, l, p, pp); 

    for (uint64_t value : pp) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    std::vector<std::vector<int>> A, B, C;
    
    try {
        lib.createMat("instructions.txt", A, B, C, p);

        // Print the matrices
        auto printMatrix = [](const std::vector<std::vector<int>>& matrix, const std::string& name) {
            std::cout << "Matrix " << name << ":" << std::endl;
            for (const auto& row : matrix) {
                for (int val : row) {
                    std::cout << val << " ";
                }
                std::cout << std::endl;
            }
        };

        printMatrix(A, "A");
        printMatrix(B, "B");
        printMatrix(C, "C");

    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}
