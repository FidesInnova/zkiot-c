// #include <cstdio>

int main() {
    int result;

    asm volatile (
        "li a0, 4\n"           // Load immediate value 4 into register a0
        "addi a0, a0, 5\n"     // Add immediate value 5 to a0 (a0 = a0 + 5)
        "mul a0, a0, 11\n"     // Multiply a0 by 11 (a0 = a0 * 11)
        "addi a0, a0, 26\n"    // Add immediate value 26 to a0 (a0 = a0 + 26)
        /*
        "mv %0, a0\n"          // Move the final result in a0 to the C++ variable `result`
        : "=r" (result)        // Output operand
        :                      // No input operands
        : "a0"                 // Clobbered register
        */
    );

    // std::cout << "Result: " << result << std::endl;

    return 0;
}