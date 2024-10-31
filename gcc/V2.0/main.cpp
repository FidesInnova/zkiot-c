#include "fidesinnova.h"

fidesinnova lib;

int main() {
  uint64_t p = 181;  // Initialize the static member 18446744069414584321 | 2013265921 | 181
  uint64_t g = 2;      // Initialize g 33 | 2
  uint64_t n = 5;      // Initialize m
  uint64_t b = 2;      // Initialize b
  uint64_t tau = 119;  // Initialize tau

  lib.setup(g, tau, p);
  // lib.commitmentGenerator("instruction.txt", g, p);
  // lib.proofGenerator("instruction.txt", g, b, p);
  // lib.verifier(g, p);
}