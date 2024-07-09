#ifndef FP_FIELD_H
#define FP_FIELD_H

#include <Arduino.h>
#include <stdint.h>

class Fp {
public:
    static uint64_t P; // Static member to hold the value of P
    uint64_t value;

    Fp();
    Fp(uint64_t v);

    static void setP(uint64_t prime); // Function to set the value of P

    // Addition
    Fp operator+(const Fp& other) const;

    // Subtraction
    Fp operator-(const Fp& other) const;

    // Multiplication
    Fp operator*(const Fp& other) const;

    // Modular exponentiation
    Fp pow(uint64_t exp) const;

    // Modular inverse using Fermat's Little Theorem
    Fp inverse() const;

    // Division
    Fp operator/(const Fp& other) const;

    // Output stream
    void print(Stream &out) const;
};

#endif
