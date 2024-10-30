#ifndef FIDESINNOVA_H
#define FIDESINNOVA_H

#include <stdint.h>
#include <Preferences.h>
#include <vector>
#include "polynomial.h"
#include <ArduinoJson.h>

using namespace std;

class fidesinnova {
public:
  void setup(int64_t g, int64_t tau, int64_t mod);
  void commitmentGenerator(String path, int64_t g, int64_t mod);
  void proofGenerator(String path, int64_t g, int64_t b, int64_t mod);
  void verifier(int64_t g, int64_t mod);

private:
};

#endif  // FIDESINNOVA_H
