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

/*
*   Arduino IDE 2.3.2
*   ESP32 by Espressif Systems
*   ESP32 Verion: 3.1.0-RC1
*   Microcontroller: ESP32-C6
*   Flash Size: 8MB
*   Partition Scheme: 8M with SPIFFS (3MB App/1.5MB SPIFFS)
*/


#include "FidesInnova.h"
#include "RegisterReader.h"
#include "SPIFFS.h"

FidesInnova lib;

void setup() {
  Serial.begin(115200);

  if (!SPIFFS.begin(true)) {
    Serial.println("SPIFFS Mount Failed");
    return;
  }

  delay(2000);
  while (!Serial) {
    delay(10);
  }

  uint64_t mod = 2147482921;  // Initialize the static member 18446744069414584321 | 2013265921 | 181 | 45151681
  uint64_t g = 2147482921;    // Initialize g 33 | 2 | 61
  uint64_t n = 5;             // Initialize m
  uint64_t b = 2;             // Initialize b
  uint64_t tau = 119;         // Initialize tau

  lib.Setup(g, tau, mod);
  lib.Commitment("instruction.txt", g, mod);
  lib.Proof("instruction.txt", g, b, mod);
  lib.Verify(g, mod);

  /*
  // Read the registers
  Registers registers = readRegisters();
  // Print the ISA value
  Serial.printf("ISA: 0x%08x (%s)\n", registers.misa, registers.readableISA.c_str());
  // Print the values of the registers
  for (int i = 0; i < 32; i++) {
    Serial.printf("X%d: 0x%08x\n", i, registers.reg_x[i]);
  }
  Serial.println("");
  // Print the values of the registers
  for (int i = 0; i < 32; i++) {
    Serial.printf("X%d: 0x%08x\n", i, registers.reg_y[i]);
  }
  */
}

void loop() {
  // put your main code here, to run repeatedly:
}
