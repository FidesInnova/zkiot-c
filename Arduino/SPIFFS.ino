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

File file;

String readFile(String fileName) {
  String Buffer = "";
  file = SPIFFS.open(fileName, "r");
  if (!file) {
    Serial.print("Failed to open ");
    Serial.println(fileName);
  } else {
    while (file.available()) {
      Buffer += file.readStringUntil('\n');
    }
  }
  file.close();
  Buffer.trim();
  return Buffer;
}


bool writeFile(String fileName, String payload) {
  bool result = false;
  if (SPIFFS.exists(fileName)) {
    SPIFFS.remove(fileName);
    delay(100);
  }
  file = SPIFFS.open(fileName, "w");
  if (!file) {
    Serial.print("Unable to create ");
    Serial.println(fileName);
    //return false if file can't be created;
  }
  int bytesWritten = file.println(payload);
  if (bytesWritten > 0) {
    Serial.print(bytesWritten);
    Serial.print(" bytes was written in ");
    Serial.println(fileName);
    result = true;
  } else {
    Serial.print("Failed to write ");
    Serial.println(fileName);
  }
  file.close();
  return result;
}


void removeFile(String fileName) {
  if (SPIFFS.exists(fileName)) SPIFFS.remove(fileName);
}
