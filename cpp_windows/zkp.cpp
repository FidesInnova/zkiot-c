#include "zkp.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include <cstring>
#include <ctype.h> 
#include <cctype>

using namespace std;


uint64_t power(uint64_t base, uint64_t exponent, uint64_t modulus) {
    uint64_t result = 1;
    base = base % modulus;
    while (exponent > 0) {
        if (exponent % 2 == 1) {
            result = (result * base) % modulus;
        }
        exponent = exponent >> 1;
        base = (base * base) % modulus;
    }
    return result;
}

void ZKP::setup(uint64_t g, uint64_t d, uint64_t l, uint64_t p, std::vector<uint64_t>& pp) {
    pp.clear();
    uint64_t pMinusOne = p - 1;
    uint64_t current_exponent = 1;
    uint64_t newPower = d % pMinusOne;
    uint64_t res = 0;

    res = power(g, newPower, p);
    pp.push_back(res);
    for (uint64_t i = 1; i < l; ++i) {
        res = power(res, newPower, p);
        pp.push_back(res);
    }
    //pp=(2,66,83,91,96,24,2,66,83)
}

// Helper function to initialize a zero matrix of given size
std::vector<std::vector<int>> initializeZeroMatrix(int size) {
    return std::vector<std::vector<int>>(size, std::vector<int>(size, 0));
}


// Helper function to trim spaces from a string
std::string trim(const std::string& str) {
    size_t first = str.find_first_not_of(' ');
    if (first == std::string::npos)
        return "";
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, last - first + 1);
}


void ZKP::createMat(const std::string& filename, std::vector<std::vector<int>>& A, std::vector<std::vector<int>>& B, std::vector<std::vector<int>>& C) {
    // Reading the instructions from the file
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Error opening file");
    }



    std::string line;
    std::vector<std::string> instructions;
    std::unordered_map<int, int> inputs; // Store which registers are inputs
    std::cout << "Reading " << filename << "\n";
    while (std::getline(file, line)) {
        std::cout << line << "\n";
        instructions.push_back(line);
    }
    std::cout << "End of the file\n";
    file.close();

    // Number of gates and inputs
    int n_g = 0;
    int n_i = 0;
    
    // Parse instructions to determine n_g and n_i
    for (const auto& instr : instructions) {
        std::stringstream ss(instr);
        std::string operation, rStr, xStr, yStr;
        
        ss >> operation >> rStr;

        if (operation == "li") {
            ss >> xStr;
            n_i++;
            int R = std::stoi(rStr.substr(1));
            int X = std::stoi(xStr);
            inputs[R] = X; // Store the input value
        } else {
            ss >> xStr >> yStr;
            n_g++;
        }
    }

    // Matrix order
    int n = n_g + n_i + 1;
    std::cout << "Matrix order: " << n << "\n";

    // Initialize matrices A, B, C
    A = initializeZeroMatrix(n);
    B = initializeZeroMatrix(n);
    C = initializeZeroMatrix(n);


    // Fill matrices based on the instructions
    int gateIndex = 0;
    for (const auto& instr : instructions) {
        std::stringstream ss(instr);
        std::string operation, rStr, xStr, yStr;
        
        ss >> operation >> rStr;
        int R = std::stoi(rStr.substr(1));

        if (operation != "li") {
            ss >> xStr >> yStr;

            // Remove commas
            xStr.pop_back();
            yStr.pop_back();

            // Trim spaces
            xStr = trim(xStr);
            yStr = trim(yStr);

            std::cout << "operation: " << operation << "\tR: " << R << "\tX: " << xStr << "\tY: " << yStr << "\n";
            // int X = std::isdigit(xStr[0]) ? std::stoi(xStr) : std::stoi(xStr.substr(1));
            // int Y = std::isdigit(yStr[0]) ? std::stoi(yStr) : std::stoi(yStr.substr(1));

            gateIndex++;
            int i = gateIndex + n_i;
            C[i][i] = 1;

            // if (operation == "add") {
            //     A[i][1] = 1;

            //     if (inputs.find(X) != inputs.end()) {
            //         B[i][1] = inputs[X]; // Set the left input value
            //     } else {
            //         B[i][1 + X] = 1; // Set the left input as a register
            //     }

            //     if (inputs.find(Y) != inputs.end()) {
            //         B[i][1] = inputs[Y]; // Set the right input value
            //     } else {
            //         B[i][1 + Y] = 1; // Set the right input as a register
            //     }
            // } else if (operation == "mul") {
            //     if (inputs.find(X) != inputs.end()) {
            //         A[i][1] = inputs[X]; // Set the left input value
            //     } else {
            //         A[i][1 + X] = 1; // Set the left input as a register
            //     }

            //     if (inputs.find(Y) != inputs.end()) {
            //         B[i][1] = inputs[Y]; // Set the right input value
            //     } else {
            //         B[i][1 + Y] = 1; // Set the right input as a register
            //     }
            // }
        }
    }
}






















//     // Fill matrices based on the instructions
//     int gateIndex = 0;
//     for (const auto& instr : instructions) {
//         std::stringstream ss(instr);
//         std::string operation;
//         char R, X, Y;
//         ss >> operation >> R >> X >> Y;

//         // parseInstruction(instr, operation, R, X, Y);
//         std::cout << "operation: " << operation << "\tR: " << R << "\tX: " << X << "\tY: " << Y << "\n";
//         if (operation == "add" || operation == "mul") {
//             gateIndex++;
//             int i = gateIndex + n_i;
//             C[i][i] = 1;
            
//             if (operation == "add") {
//                 A[i][1] = 1;

//                 // if (inputs.find(X) != inputs.end()) {
//                 //     B[i][1] = inputs[X]; // Set the left input value
//                 // } else {
//                 //     B[i][1 + X] = 1; // Set the left input as a register
//                 // }

//                 // if (inputs.find(Y) != inputs.end()) {
//                 //     B[i][1] = inputs[Y]; // Set the right input value
//                 // } else {
//                 //     B[i][1 + Y] = 1; // Set the right input as a register
//                 // }
//             } else if (operation == "mul") {
//                 // if (inputs.find(X) != inputs.end()) {
//                 //     A[i][1] = inputs[X]; // Set the left input value
//                 // } else {
//                 //     A[i][1 + X] = 1; // Set the left input as a register
//                 // }

//                 // if (inputs.find(Y) != inputs.end()) {
//                 //     B[i][1] = inputs[Y]; // Set the right input value
//                 // } else {
//                 //     B[i][1 + Y] = 1; // Set the right input as a register
//                 // }
//             }
//         }
//     }
//     std::cout << "gateIndex: " << gateIndex << "\n";
// }
