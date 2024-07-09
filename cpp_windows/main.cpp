#include <iostream>
#include "zkp.h"

int main() {
    uint64_t p = 181;
    uint64_t g = 2;
    uint64_t d = 111213119;
    uint64_t l = 8;
    std::vector<uint64_t> pp;

    ZKP lib;
    lib.setup(g, d, l, p, pp); 

    for (uint64_t value : pp) {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    return 0;
}
