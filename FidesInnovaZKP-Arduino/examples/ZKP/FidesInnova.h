#ifndef FIDESINNOVA_H
#define FIDESINNOVA_H

#include <stdint.h>
#include <Arduino.h>
#include <FS.h>
#include <SPIFFS.h>
#include <Preferences.h>
#include <vector>
#include "polynomial.h"

using namespace std;

class FidesInnova {
public:
  // void Setup(int64_t g, int64_t d, int64_t l, int64_t mod);
  void Commitment(String path, int64_t g, int64_t b, int64_t mod);
  void Proof();
  void Verify();

private:
};

#endif  // FIDESINNOVA_H
