#include <iostream>
#include <fstream>
#include "json.hpp"

using json = nlohmann::json;

int main() {
    std::cout << "Current directory: " << std::ifstream::current_path() << std::endl;

    // Open the JSON file
    std::ifstream file("setup3.json");
    if (!file.is_open()) {
        std::cerr << "Error opening file!" << std::endl;
        return 1;
    }

    // Parse the JSON content
    json j;
    file >> j;

    // Print the content of the JSON
    std::cout << "JSON Content: " << std::endl;
    std::cout << j.dump(4) << std::endl; // Pretty print with an indentation of 4 spaces

    return 0;
}
