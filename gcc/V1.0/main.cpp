#include <iostream>
#include <vector>
#include "zkp.h"

int main() {
    int64_t p = 181;// 18446744069414584321;// 2013265921; // 181;
    int64_t g = 2;
    int64_t d = 111213119;
    int64_t l = 8;
    int64_t n = 5;
    std::vector<int64_t> pp;

    ZKP lib;
    lib.setup(g, d, l, p, pp);

    std::cout << "pp: ";
    for (int64_t value : pp) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    std::vector<std::vector<int64_t>> A, B, C;
    
    try {
        lib.createMat("instructions.txt", A, B, C, p);

        // Print the matrices
        auto printMatrix = [](const std::vector<std::vector<int64_t>>& matrix, const std::string& name) {
             std::cout << "Matrix " << name << ":" << std::endl;
             for (const auto& row : matrix) {
                 for (int64_t val : row) {
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

    lib.AHPverify(p);
    
    return 0;
}
