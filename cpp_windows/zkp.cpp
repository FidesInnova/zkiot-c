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

// Helper function to Remove commas from a string
std::string removeCommas(const std::string& str) {
    size_t first = str.find_first_not_of(',');
    if (first == std::string::npos)
        return "";
    size_t last = str.find_last_not_of(',');
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
            xStr = trim(xStr);
            yStr = trim(yStr);
            xStr = removeCommas(xStr);
            yStr = removeCommas(yStr);
            n_g++;
        }
        std::cout << "operation: " << operation << "\tX: " << xStr << "\tY: " << yStr << "\n";
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
        int64_t z[n];
        z[0] = 1;
    for (int i = 0; i < n-1; i++) {
        std::stringstream ss(instructions[i]);
        std::string operation, rStr, xStr, yStr;
        ss >> operation >> rStr;
        int R = std::stoi(rStr.substr(1));

        if (operation == "li") {
            ss >> xStr;
        }
        else {
            ss >> xStr >> yStr;
        }


        // Remove commas
        xStr = removeCommas(xStr);
        yStr = removeCommas(yStr);

        // Trim spaces
        xStr = trim(xStr);
        yStr = trim(yStr);

        std::cout << "operation: " << operation << "\tX: " << xStr << "\tY: " << yStr << "\n";

        // int X = std::isdigit(xStr[0]) ? std::stoi(xStr) : std::stoi(xStr.substr(1));
        // int Y = std::isdigit(yStr[0]) ? std::stoi(yStr) : std::stoi(yStr.substr(1));
        int X,Y;
        if (operation == "li") {
            X = std::stoi(xStr);
            z[i+1] = X;
        } else {
            ss >> xStr >> yStr;
            
            gateIndex++;
            int i = gateIndex + n_i;
            C[i][i] = 1;

            if (operation == "add") {
                A[n_i+i-1][1-1] = 1;

                if (std::isdigit(xStr[0])) {
                    X = std::stoi(xStr);
                    z[i+1] = z[i] + X;
                    B[n_i+i-1][0] = X;
                    B[n_i+i-1][i-1] = 1;
                } else if(std::isdigit(yStr[0])){
                    Y = std::stoi(yStr);
                    z[i+1] = z[i-1] + Y;
                    B[n_i+i-1][0] = Y;
                    B[n_i+i-1][i-1] = 1;
                }

            } else if (operation == "mul") {
                if (std::isdigit(xStr[0])) {
                    X = std::stoi(xStr);
                    z[i+1] = z[i] * X;
                    A[n_i+i-1][i-1] = X;
                    B[n_i+i-1][i-1] = 1;
                } else if (std::isdigit(yStr[0])) {
                    Y = std::stoi(yStr);
                    z[i+1] = z[i] * Y;
                    A[n_i+i-1][i-1] = 1;
                    B[n_i+i-1][0] = Y;
                }
            }
        }
    }
}

