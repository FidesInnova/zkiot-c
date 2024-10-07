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
    //return;
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
