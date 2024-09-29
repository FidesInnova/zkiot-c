#ifndef RegisterReader_h
#define RegisterReader_h

#include <Arduino.h>
#include <cstdint>

// Define a structure to hold registers, PC, ISA, and readable ISA
struct Registers {
  uint32_t pc;
  uint32_t reg_x[32];
  uint32_t reg_y[32];
  uint32_t misa;
  String readableISA;
};

Registers readRegisters();

#endif
