#include <stdint.h>
#include <fstream>
#include "json.hpp"
using ordered_json = nlohmann::ordered_json;
#include <regex>
#include <iostream>

using namespace std;

void setup(int64_t tau) {
    vector<int64_t> ck;
    int64_t vk;

    nlohmann::json classJsonData;
    // Open the JSON file for reading (class.json)
    std::ifstream classFile("class.json");
    if (!classFile.is_open()) {
        cerr << "Could not open class.json!" << endl;
        return; // Exit if the file cannot be opened
    }

    // Parse the JSON data from class.json
    classFile >> classJsonData;
    classFile.close();

    // Loop to create multiple setup files
    for (int class_value = 1; class_value <= 16; class_value++) {
        string Class = to_string(class_value); // Convert integer to string class
        // Get n_g for current file
        if (classJsonData.contains(Class)) { // Check if the class exists
            int64_t n_g = classJsonData[Class]["n_g"].get<int64_t>();
            int64_t n_i = classJsonData[Class]["n_i"].get<int64_t>();
            int64_t n   = classJsonData[Class]["n"].get<int64_t>();
            int64_t m   = classJsonData[Class]["m"].get<int64_t>();
            int64_t p   = classJsonData[Class]["p"].get<int64_t>();
            int64_t g   = classJsonData[Class]["g"].get<int64_t>();

            cout << "Class: " << class_value << endl;
            cout << "  n_g: " << n_g << endl;
            cout << "  n_i: " << n_i << endl;
            cout << "  n: " << n << endl;
            cout << "  m: " << m << endl;
            cout << "  p: " << p << endl;
            cout << "  g: " << g << endl;

            // Validate p and g
            if (p <= 1) {
                cout << "Invalid p value for Class " << class_value << ": " << p << endl;
            } else {
                int64_t d_AHP = 12 * n_g;

                // Calculate new power for g based on d_AHP
                int64_t pMinusOne = p - 1;
                tau %= pMinusOne;

                // Push values into ck
                ck.clear(); // Clear previous values

                for (int64_t i = 0; i < d_AHP; i++) {
                    ck.push_back(g);
                    g = (g * tau) %  p;
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
                setupJson["Class"] = class_value;
                setupJson["ck"] = ck;
                setupJson["vk"] = vk;

                // Serialize JSON object to a string
                string setupString = setupJson.dump(4); // Pretty print with 4-space indentation

                // Write JSON object to a file
                std::ofstream setupFile("data/setup" + to_string(class_value) + ".json");
                if (setupFile.is_open()) {
                    setupFile << setupString;
                    setupFile.close();
                    cout << "JSON data has been written to setup" << class_value << ".json\n";
                } else {
                    cerr << "Error opening file for writing setup" << class_value << ".json\n";
                }
            }
        } else {
            cout << "Class " << class_value << " not found in JSON.\n";
        }
    }
}

int main() {
    int64_t tau = 119; // Initialize tau
    setup(tau); // Call setup with tau
}
