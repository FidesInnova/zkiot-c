#include "FpField.h"

uint64_t Fp::P = 1; // Initialize the static member

Fp::Fp() : value(0) {}

Fp::Fp(uint64_t v) : value(v % P) {}

void Fp::setP(uint64_t prime) {
    P = prime;
}

// Addition
Fp Fp::operator+(const Fp& other) const {
    return Fp((value + other.value) % P);
}

// Subtraction
Fp Fp::operator-(const Fp& other) const {
    return Fp((value + P - other.value) % P);
}

// Multiplication
Fp Fp::operator*(const Fp& other) const {
    return Fp((value * other.value) % P);
}

// Modular exponentiation
Fp Fp::pow(uint64_t exp) const {
    uint64_t result = 1;
    uint64_t base = value;
    while (exp > 0) {
        if (exp % 2 == 1) {
            result = (result * base) % P;
        }
        exp = exp >> 1;
        base = (base * base) % P;
    }
    return Fp(result);
}

// Modular inverse using Fermat's Little Theorem
Fp Fp::inverse() const {
    return pow(P - 2);
}

// Division
Fp Fp::operator/(const Fp& other) const {
    return *this * other.inverse();
}

// Output stream
void Fp::print(Stream &out) const {
    out.print(value);
}
