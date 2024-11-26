#include "fidesinnova.h"

int main() {
    int result;

    asm volatile (
        // "li s1, 4\n"
        // "li s1, 5\n"
        // "li s1, 26\n"
        // "add s1, s1, s1\n"
        // "addi s1, s1, 11\n"
        // "add s1, s1, s1\n"
        // "add s1, s1, s1\n"

        // "li s1, 1\n"
        // "add s1, s1, s1\n"
        // "addi s1, s1, 11\n"
        // "add s1, s1, s1\n"
        // "add s1, s1, s1\n"
        // "add s1, s1, s1\n"
        // "add s1, s1, s1\n"
        // "add s1, s1, s1\n"
        // "add s1, s1, s1\n"

        
        "li s1, 1\n"
        "li s2, 2\n"
        "li s3, 3\n"
        "li s4, 4\n"
        "li s5, 5\n"
        "mul s2, s2, s3\n"
        "addi s1, s3, 11\n"
        "mul s2, s4, s4\n"
        "mul s3, s3, s4\n"
        "add s1, s2, s3\n"
        "mul s3, s4, s4\n"
        "mul s2, s1, s1\n"
        "mul s1, s2, s2\n"
        "mul s2, s2, s3\n"
        "addi s1, s3, 11\n"
        "mul s2, s4, s4\n"
        "mul s3, s3, s4\n"
        "add s1, s2, s3\n"
        "mul s3, s4, s4\n"
        "mul s2, s1, s1\n"
        "mul s1, s2, s2\n"
    );
    return 0;
}