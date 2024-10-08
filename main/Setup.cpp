#include "FidesInnova.h"
#include <stdint.h>


void FidesInnova::Setup(int64_t g, int64_t tau, int64_t mod) {
  vector<int64_t> ck;

  /*
  // Calculate each expression
  int64_t exp1 = m;
  int64_t exp2 = n_g - n_i + b;
  int64_t exp3 = n + b;
  int64_t exp4 = n + 2 * b - 1;
  int64_t exp5 = 2 * n + b - 1;
  int64_t exp6 = n + b - 1;
  int64_t exp7 = n - 1;
  int64_t exp8 = m - 1;
  int64_t exp9 = 6 * m - 6;

  // Find the maximum value
  int64_t d_AHP = max(exp1, max(exp2, max(exp3, max(exp4, max(exp5, max(exp6, max(exp7, max(exp8, exp9))))))));
  */
  // Use a predefined max
  int64_t d_AHP = 100;

  int64_t modMinusOne = mod - 1;
  int64_t current_exponent = 1;
  int64_t newPower = d_AHP % modMinusOne;
  int64_t tmp = 0;

  tau %= mod - 1;
  tmp = g;
  for (int64_t i = 0; i < d_AHP; i++) {
    ck.push_back(tmp);
    tmp = (tmp * tau) % mod;
  }

  Serial.print("ck = {");
  for (int64_t i = 0; i < ck.size(); i++) {
    Serial.print(String(ck[i]));
    if (ck.size() - i > 1) {
      Serial.print(", ");
    }
  }
  Serial.println("}");
  int64_t vk = ck[1];
  Serial.print("vk = ");
  Serial.println(String(vk));

  DynamicJsonDocument doc(2048);
  JsonArray jsonArray = doc.createNestedArray("ck");
  for (int64_t value : ck) {
    jsonArray.add(value);  // Add each element to the JSON array
  }
  jsonArray = doc.createNestedArray("vk");
  jsonArray.add(vk);
  // Serialize the JSON document to a string for testing
  String output;
  serializeJson(doc, output);

  // Store setup(ck, vk) in setup.json
  removeFile("/setup.json");
  writeFile("/setup.json", output);
}
