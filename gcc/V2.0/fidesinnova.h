#ifndef FIDESINNOVA_H
#define FIDESINNOVA_H

#include <stdint.h>
#include "polynomial.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>

using namespace std;

class fidesinnova {
public:
  void setup(int64_t tau);
  // void commitmentGenerator(String path, int64_t g, int64_t p);
  // void proofGenerator(String path, int64_t g, int64_t b, int64_t p);
  // void verifier(int64_t g, int64_t p);

private:
};

#endif  // FIDESINNOVA_H
