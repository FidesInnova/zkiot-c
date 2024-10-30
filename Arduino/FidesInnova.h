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
  void Commitment(String path, int64_t g, int64_t mod);
  void Proof(String path, int64_t g, int64_t b, int64_t mod);
  void Verify(int64_t g, int64_t mod);

private:
};

#endif  // FIDESINNOVA_H
