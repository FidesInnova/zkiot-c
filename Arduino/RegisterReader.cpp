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

#include "RegisterReader.h"

Registers readRegisters() {
  Registers result;

  // Inline assembly to get the value of the PC

  // __asm__ volatile ("jalr %0, 0" : "=r"(result.pc));    // Get current PC
  // __asm__ volatile (
  //     "auipc %0, 0\n"        // Load upper immediate to PC
  //     "addi  %0, %0, 8\n"    // Add the immediate value to get the address of the next instruction
  //     : "=r" (result.pc)
  // );
  // __asm__ volatile (
  //     "auipc t0, 0\n"        // Load upper immediate to PC
  //     "addi  t0, t0, 8\n"    // Calculate the address of the next instruction
  //     "jalr  x0, t0, 0\n"    // Jump to that address, x0 means no link (no return address stored)
  //     "addi  %0, t0, 0\n"    // Store the address in pc (just to demonstrate)
  //     : "=r" (result.pc)
  //     :
  //     : "t0"
  // );

  // Inline assembly to get the value of the misa register
  __asm__ volatile("csrr %0, misa"
                   : "=r"(result.misa));  // Get value of misa register

  // Interpret the ISA from misa register
  result.readableISA = "";
  if (result.misa & (1 << ('I' - 'A'))) result.readableISA += "I";
  if (result.misa & (1 << ('M' - 'A'))) result.readableISA += "M";
  if (result.misa & (1 << ('A' - 'A'))) result.readableISA += "A";
  if (result.misa & (1 << ('F' - 'A'))) result.readableISA += "F";
  if (result.misa & (1 << ('D' - 'A'))) result.readableISA += "D";
  if (result.misa & (1 << ('C' - 'A'))) result.readableISA += "C";
  // Add other ISA extensions as needed



  __asm__ volatile("li x18, 10\n");
  __asm__ volatile("li x19, 20\n");


  // Using multiple inline assembly blocks to get the values of the registers
  __asm__ volatile(
    "mv %0, x0\n"
    "mv %1, x1\n"
    "mv %2, x2\n"
    "mv %3, x3\n"
    "mv %4, x4\n"
    "mv %5, x5\n"
    "mv %6, x6\n"
    "mv %7, x7\n"
    : "=r"(result.reg_x[0]), "=r"(result.reg_x[1]), "=r"(result.reg_x[2]), "=r"(result.reg_x[3]), "=r"(result.reg_x[4]), "=r"(result.reg_x[5]), "=r"(result.reg_x[6]), "=r"(result.reg_x[7]));

  __asm__ volatile(
    "mv %0, x8\n"
    "mv %1, x9\n"
    "mv %2, x10\n"
    "mv %3, x11\n"
    "mv %4, x12\n"
    "mv %5, x13\n"
    "mv %6, x14\n"
    "mv %7, x15\n"
    : "=r"(result.reg_x[8]), "=r"(result.reg_x[9]), "=r"(result.reg_x[10]), "=r"(result.reg_x[11]), "=r"(result.reg_x[12]), "=r"(result.reg_x[13]), "=r"(result.reg_x[14]), "=r"(result.reg_x[15]));

  __asm__ volatile(
    "mv %0, x16\n"
    "mv %1, x17\n"
    "mv %2, x18\n"
    "mv %3, x19\n"
    "mv %4, x20\n"
    "mv %5, x21\n"
    "mv %6, x22\n"
    "mv %7, x23\n"
    : "=r"(result.reg_x[16]), "=r"(result.reg_x[17]), "=r"(result.reg_x[18]), "=r"(result.reg_x[19]), "=r"(result.reg_x[20]), "=r"(result.reg_x[21]), "=r"(result.reg_x[22]), "=r"(result.reg_x[23]));

  __asm__ volatile(
    "mv %0, x24\n"
    "mv %1, x25\n"
    "mv %2, x26\n"
    "mv %3, x27\n"
    "mv %4, x28\n"
    "mv %5, x29\n"
    "mv %6, x30\n"
    "mv %7, x31\n"
    : "=r"(result.reg_x[24]), "=r"(result.reg_x[25]), "=r"(result.reg_x[26]), "=r"(result.reg_x[27]), "=r"(result.reg_x[28]), "=r"(result.reg_x[29]), "=r"(result.reg_x[30]), "=r"(result.reg_x[31]));

  __asm__ volatile("add x20, x18, x19\n");




  // Using multiple inline assembly blocks to get the values of the registers
  __asm__ volatile(
    "mv %0, x0\n"
    "mv %1, x1\n"
    "mv %2, x2\n"
    "mv %3, x3\n"
    "mv %4, x4\n"
    "mv %5, x5\n"
    "mv %6, x6\n"
    "mv %7, x7\n"
    : "=r"(result.reg_y[0]), "=r"(result.reg_y[1]), "=r"(result.reg_y[2]), "=r"(result.reg_y[3]), "=r"(result.reg_y[4]), "=r"(result.reg_y[5]), "=r"(result.reg_y[6]), "=r"(result.reg_y[7]));

  __asm__ volatile(
    "mv %0, x8\n"
    "mv %1, x9\n"
    "mv %2, x10\n"
    "mv %3, x11\n"
    "mv %4, x12\n"
    "mv %5, x13\n"
    "mv %6, x14\n"
    "mv %7, x15\n"
    : "=r"(result.reg_y[8]), "=r"(result.reg_y[9]), "=r"(result.reg_y[10]), "=r"(result.reg_y[11]), "=r"(result.reg_y[12]), "=r"(result.reg_y[13]), "=r"(result.reg_y[14]), "=r"(result.reg_y[15]));

  __asm__ volatile(
    "mv %0, x16\n"
    "mv %1, x17\n"
    "mv %2, x18\n"
    "mv %3, x19\n"
    "mv %4, x20\n"
    "mv %5, x21\n"
    "mv %6, x22\n"
    "mv %7, x23\n"
    : "=r"(result.reg_y[16]), "=r"(result.reg_y[17]), "=r"(result.reg_y[18]), "=r"(result.reg_y[19]), "=r"(result.reg_y[20]), "=r"(result.reg_y[21]), "=r"(result.reg_y[22]), "=r"(result.reg_y[23]));

  __asm__ volatile(
    "mv %0, x24\n"
    "mv %1, x25\n"
    "mv %2, x26\n"
    "mv %3, x27\n"
    "mv %4, x28\n"
    "mv %5, x29\n"
    "mv %6, x30\n"
    "mv %7, x31\n"
    : "=r"(result.reg_y[24]), "=r"(result.reg_y[25]), "=r"(result.reg_y[26]), "=r"(result.reg_y[27]), "=r"(result.reg_y[28]), "=r"(result.reg_y[29]), "=r"(result.reg_y[30]), "=r"(result.reg_y[31]));

  return result;
}
