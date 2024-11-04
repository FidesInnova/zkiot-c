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

  // Define an array for variable n_g values for the 12 setups
  vector<int> n_g_values = {1,2,4,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536};

  // Loop to create multiple setup files
  for (int class_value = 0; class_value < n_g_values.size(); ++class_value) {
    // Get n_g for current file
    int n_g = n_g_values[class_value];
    int64_t d_AHP = 12 * n_g; // Calculate d_AHP based on n_g

    // Calculate new power for g based on d_AHP
    int64_t pMinusOne = p - 1;
    tau %= pMinusOne;

    // Push values into ck
    ck.clear(); // Clear previous values
    for (int64_t i = 0; i < d_AHP; i++) {
        ck.push_back(g);
        g = (g * tau) % p;
    }

    // Output ck for verification
    cout << "ck = {";
    for (int64_t i = 0; i < ck.size(); i++) {
        cout << ck[i];
        if (ck.size() - i > 1) {
            cout << ", ";
        }
    }
    cout << "}" << endl;

    // Retrieve verifying key
    if (ck.size() > 1) {
        vk = ck[1];
        cout << "vk = " << vk << endl;
    } else {
        cout << "Error: ck does not have enough elements." << endl;
        vk = 0; // Set vk to 0 if ck is insufficient
    }

    // Create JSON object
    ordered_json setupJson;
    setupJson["Class"] = class_value + 1; // Class number starting from 1
    setupJson["ck"] = ck;
    setupJson["vk"] = vk;

    // Serialize JSON object to a string
    std::string setupString = setupJson.dump(4); // Pretty print with 4-space indentation

    // Write JSON object to a file
    std::ofstream setupFile("data/setup" + std::to_string(class_value + 1) + ".json");
    if (setupFile.is_open()) {
        setupFile << setupString;
        setupFile.close();
        cout << "JSON data has been written to setup" << class_value + 1 << ".json\n";
    } else {
        cerr << "Error opening file for writing setup" << class_value + 1 << ".json\n";
    }
  }
}
