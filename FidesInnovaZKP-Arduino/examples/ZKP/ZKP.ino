#include "FidesInnova.h"
FidesInnova lib;

void setup() {
  Serial.begin(115200);
  uint64_t mod = 181;      // Initialize the static member 18446744069414584321 | 2013265921 | 181
  uint64_t g = 2;          // Initialize g 33 | 2
  uint64_t d = 111213119;  // Initialize d
  uint64_t l = 8;          // Initialize l
  uint64_t n = 5;          // Initialize m
  uint64_t b = 2;          // Initialize b

  lib.Setup(g, d, l, mod);
  lib.Commitment("instruction.txt", g, b, mod);
}

void loop() {
  // put your main code here, to run repeatedly:
}
