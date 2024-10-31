#include "fidesinnova.h"
#include <stdint.h>

#include "json.hpp"
using ordered_json = nlohmann::ordered_json;

#include <regex>

#include <iostream>
using namespace std;

void fidesinnova::setup(int64_t g, int64_t tau, int64_t p) {
  vector<int64_t> ck;
  int64_t vk;

  // TheArch->TheSol: Test and uncomment the following section.
  /*
  // //I add this to ensure the value of p is valid, or it will show overflow error
  // if (p > 1) { 
  //   tau %= p - 1;
  // } else {
  //   cout << "Invalid pulus value" << endl;
  //   return;
  // } */

  // TheArch->TheSol: Test and uncomment the following section.
  // b<n_g, the number of elements in c_k < 12*n_g, if verificaion returns true and double-check confirms that the number of non-zero elements in matrix B is 2.
  // TheArch->TheSol: MAke sure b <n_g in the proof generation phase.
  // TheArch->Kath: Generate 12 setup.json with the maximum of c_k elements of 12*n_g for each class. 
  // Class defenition is here: https://app.gitbook.com/o/viOQScpcYEhrVQAFsW8c/s/6jhloRdZDsdxHkgxETHg/zero-knowledge-proof-zkp-scheme/1-setup-phase
  // setupX.json format is here:
  // {
  //  "Class":  32-bit Integer,
  //  "ck": 64-bit Integer Array,
  //  "vk": 64-bit Integer
  // }
  
  // Find the maximum value
  // int64_t d_AHP = 12 * n_g;
  
  // Use a predefined max
  int64_t d_AHP = 100;

  int64_t pMinusOne = p - 1;
  int64_t current_exponent = 1;
  int64_t newPower = d_AHP % pMinusOne;

  tau %= p - 1;

  //push into ck
  for (int64_t i = 0; i < d_AHP; i++) {
    ck.push_back(g);
    g = (g * tau) % p;

  }

  cout << "ck = {";
  for (int64_t i = 0; i < ck.size(); i++) {
    cout << ck[i];
    if (ck.size() - i > 1) {
      cout << ", ";
    }
  }
  cout << endl;

  //retrive verifying key
  if (ck.size() > 1) {
    vk = ck[1];
    cout << "vk = " << vk << endl;
  } else {
    cout << "Error: ck does not have enough elements." << endl;
  }

  // TODO: Add class numbert to the setup.json file

  /*
  //it might be more clear to use different name to store JsonArray between "ck" and "vk", 
  // as the reference to the first array will be overwritten when the second array is created.
  // Create the second nested array under "vk"
  // JsonArray vkArray = doc.createNestedArray("vk");
  // vkArray.add(vk);  // Add the "vk" value to the "vk" JSON array */

  ordered_json setupJson;
  setupJson.clear();
  setupJson["ck"] = ck;
  setupJson["vk"] = vk;

  // Serialize JSON object to a string
  std::string setupString = setupJson.dump();

  // Use regex to format arrays correctly
  // Write JSON object to a file
  std::ofstream setupFile("data/setupX.json");
  if (setupFile.is_open()) {
      setupFile << setupString;
      setupFile.close();
      cout << "JSON data has been written to setupX.json\n";
  } else {
      cerr << "Error opening file for writing\n";
  }

  // TheArch->Kath: Generate 18 setup.json with the maximum of c_k elements of 12*n_g for each class. 
  // setup1.json, setup2.json, setup3.json, ..., setup18.json
  // Store setup(class, ck, vk) in setupX.json
}
