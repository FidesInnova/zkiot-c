// Copyright 2024 FidesInnova.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

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
