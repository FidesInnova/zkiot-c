#include "fidesinnova.h"
#include <stdint.h>


void fidesinnova::setup(int64_t g, int64_t tau, int64_t mod) {
  vector<int64_t> ck;

  /*
  // //I add this to ensure the value of mod is valid, or it will show overflow error
  // if (mod > 1) { 
  //   tau %= mod - 1;
  // } else {
  //   Serial.println("Invalid modulus value");
  //   return;
  // } */


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

  /*
  //can directly assign tmp to g
  // int64_t tmp = g;*/

  int64_t tmp = 0;

  tau %= mod - 1;
  tmp = g;

  //push into ck
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

  //retrive verifying key
  Serial.println("}");
  int64_t vk = ck[1];
  Serial.print("vk = ");
  Serial.println(String(vk));

  // //I add this to check ck size is bigger than 1, as it will avoid overflow problem
  // if (ck.size() > 1) {
  //   int64_t vk = ck[1];
  //   Serial.print("vk = ");
  //   Serial.println(String(vk));
  // } else {
  //   Serial.println("Error: ck does not have enough elements.");
  // }

  //Memory consideration
  DynamicJsonDocument doc(2048);

  // TODO: Add class numbert to the setup.json file



  JsonArray jsonArray = doc.createNestedArray("ck");
  for (int64_t value : ck) {
    jsonArray.add(value);  // Add each element to the JSON array
  }

  /*
  //it might be more clear to use different name to store JsonArray between "ck" and "vk", 
  // as the reference to the first array will be overwritten when the second array is created.
  // Create the second nested array under "vk"
  // JsonArray vkArray = doc.createNestedArray("vk");
  // vkArray.add(vk);  // Add the "vk" value to the "vk" JSON array */

  jsonArray = doc.createNestedArray("vk");
  jsonArray.add(vk);
  // Serialize the JSON document to a string for testing
  String output;
  serializeJson(doc, output);

  // Store setup(ck, vk) in setup.json
  removeFile("/setup.json");
  writeFile("/setup.json", output);

}
