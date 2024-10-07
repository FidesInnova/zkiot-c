#ifndef FIDESINNOVA_H
#define FIDESINNOVA_H

#include <stdint.h>
#include <Arduino.h>
#include <FS.h>
#include <SPIFFS.h>
#include <Preferences.h>
#include <vector>
#include "polynomial.h"
#include <ArduinoJson.h>

String readFile(String fileName);
bool writeFile(String fileName, String payload);
void removeFile(String fileName);

using namespace std;

class FidesInnova {
public:
  void Setup(int64_t g, int64_t tau, int64_t mod);
  void Commitment(String path, int64_t g, int64_t b, int64_t mod);
  void Proof(String path, int64_t g, int64_t b, int64_t mod);
  void Verify(int64_t mod);

private:
};

#endif  // FIDESINNOVA_H
