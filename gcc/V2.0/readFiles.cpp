#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "json.hpp"

// Function to read JSON config file and parse lines to read from assembly file
std::pair<int, int> parseDeviceConfig(const std::string &configFile, nlohmann::json &config) {
    std::ifstream configFileStream(configFile, std::ifstream::binary);
    if (!configFileStream.is_open()) {
        std::cerr << "Error opening config file: " << configFile << std::endl;
        exit(EXIT_FAILURE);
    }

    configFileStream >> config;
    configFileStream.close();

    std::vector<int> linesToRead;

    int startLine = config["Lines"][0].get<int32_t>();
    int endLine = config["Lines"][1].get<int32_t>();

    return {startLine, endLine};
}

// Function to read specified lines from assembly file
std::vector<std::string> readAssemblyLines(const std::string &assemblyFile, int startLine, int endLine) {
    std::ifstream assemblyFileStream(assemblyFile);
    if (!assemblyFileStream.is_open()) {
        std::cerr << "Error opening assembly file: " << assemblyFile << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<std::string> selectedLines;
    std::string line;
    int currentLineNumber = 1;

    while (std::getline(assemblyFileStream, line)) {
        if (currentLineNumber >= startLine && currentLineNumber <= endLine) {
            selectedLines.push_back(line);
        }
        ++currentLineNumber;
    }

    assemblyFileStream.close();
    return selectedLines;
}

// Function to modify assembly file content and save to new file
void modifyAndSaveAssembly(const std::string &assemblyFile, const std::string &newAssemblyFile, int startLine, int endLine) {
    std::ifstream assemblyFileStream(assemblyFile);
    if (!assemblyFileStream.is_open()) {
        std::cerr << "Error opening assembly file: " << assemblyFile << std::endl;
        exit(EXIT_FAILURE);
    }

    std::ofstream newAssemblyFileStream(newAssemblyFile);
    if (!newAssemblyFileStream.is_open()) {
        std::cerr << "Error creating new assembly file: " << newAssemblyFile << std::endl;
        exit(EXIT_FAILURE);
    }

    std::string line;
    int currentLineNumber = 1;
    int index = 0;

    while (std::getline(assemblyFileStream, line)) {
        if (currentLineNumber == 1) {
            newAssemblyFileStream << ".data               # Data section\n";
            newAssemblyFileStream << "X: .word 0         # Reserve space for X\n";
            for(int i = 0; i < endLine - startLine; i++) {
                std::string w_data = "W" + std::to_string(i) + ": .word 0        # Reserve space for W" + std::to_string(i) + "\n";
                newAssemblyFileStream << w_data;
            }
            newAssemblyFileStream << "Y: .word 0         # Reserve space for Y\n";
        }
        // Insert variables before the specified lines
        if (currentLineNumber == startLine) {
            newAssemblyFileStream << "sw a0, X\n";
        }
        else if (currentLineNumber == endLine + 1){
            newAssemblyFileStream << "sw a0, Y\n";
        }
        else if(currentLineNumber > startLine && currentLineNumber <= endLine){
            std::string w = "sw a0, W" + std::to_string(currentLineNumber - (startLine + 1)) + "\n";
            newAssemblyFileStream << w;
        }

        newAssemblyFileStream << line << std::endl;
        ++currentLineNumber;
    }

    assemblyFileStream.close();
    newAssemblyFileStream.close();
}

int main() {
    std::string configFile, assemblyFile, newAssemblyFile;

    // Input filenames
    std::cout << "Enter the device config file name: ";
    std::cin >> configFile;
    std::cout << "Enter the program assembly file name: ";
    std::cin >> assemblyFile;
    std::cout << "Enter the output file name for modified assembly: ";
    std::cin >> newAssemblyFile;

    nlohmann::json config;
    auto [startLine, endLine] = parseDeviceConfig(configFile, config);

    modifyAndSaveAssembly(assemblyFile, newAssemblyFile, startLine, endLine);

    std::cout << "Modified assembly file saved as: " << newAssemblyFile << std::endl;

    return 0;
}
