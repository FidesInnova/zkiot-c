#include "zkp.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include <cstring>
#include <ctype.h> 
#include <cctype>

#include <tuple>

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

// Function to calculate the modular inverse using Extended Euclidean Algorithm
uint64_t mod_inverse(uint64_t a, uint64_t mod) {
    int64_t t = 0, newt = 1;
    int64_t r = mod, newr = a;

    while (newr != 0) {
        int64_t quotient = r / newr;
        t = t - quotient * newt;
        std::swap(t, newt);
        r = r - quotient * newr;
        std::swap(r, newr);
    }

    if (r > 1) return 0; // a is not invertible
    if (t < 0) t += mod;
    return t;
}

// Function to calculate the Lagrange basis polynomial
uint64_t lagrange_basis(const std::vector<uint64_t>& x, int i, uint64_t xi, uint64_t mod) {
    uint64_t result = 1;
    int n = x.size();
    for (int j = 0; j < n; j++) {
        if (j != i) {
            uint64_t numerator = (xi + x[j]) % mod;
            uint64_t denominator = (x[i] + x[j]) % mod;
            uint64_t inv_denominator = mod_inverse(denominator, mod);
            result = (result * numerator % mod) * inv_denominator % mod;
        }
    }
    return result;
}

// Function to interpolate using Lagrange polynomial
std::tuple<uint64_t, uint64_t, uint64_t> lagrange_interpolation(const std::vector<uint64_t>& x, const std::vector<uint64_t>& y, uint64_t mod) {
    int n = x.size();
    uint64_t c2 = 0, c1 = 0, c0 = 0;

    for (int i = 0; i < n; i++) {
        uint64_t L_i2 = lagrange_basis(x, i, 2, mod); // coefficient of x^2
        uint64_t L_i1 = lagrange_basis(x, i, 1, mod); // coefficient of x
        uint64_t L_i0 = lagrange_basis(x, i, 0, mod); // constant term

        c2 = (c2 + y[i] * L_i2 % mod) % mod;
        c1 = (c1 + y[i] * L_i1 % mod) % mod;
        c0 = (c0 + y[i] * L_i0 % mod) % mod;
    }

    return {c2, c1, c0};
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


void ZKP::createMat(const std::string& filename, std::vector<std::vector<int>>& A, std::vector<std::vector<int>>& B, std::vector<std::vector<int>>& C, uint64_t p) {
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
        // std::cout << line << "\n";
        instructions.push_back(line);
    }
    // std::cout << "End of the file\n";
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
            inputs[R] = X;
        } else {
            ss >> xStr >> yStr;
            xStr = trim(xStr);
            yStr = trim(yStr);
            xStr = removeCommas(xStr);
            yStr = removeCommas(yStr);
            n_g++;
        }
        // std::cout << "operation: " << operation << "\tX: " << xStr << "\tY: " << yStr << "\n";
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
    int64_t z[n+1];
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

        int X,Y;
        if (operation == "li") {
            X = std::stoi(xStr);
            z[i+1] = X % p;
        } else {
            ss >> xStr >> yStr;
            
            gateIndex++;
            int newI = gateIndex + n_i;
            C[newI][newI] = 1;

            if (operation == "add") {
                A[n_i+newI-1][1-1] = 1;

                if (std::isdigit(xStr[0])) {
                    X = std::stoi(xStr);
                    z[i+1] = (z[i] + X) % p;
                    B[n_i+newI-1][0] = X;
                    B[n_i+newI-1][newI-1] = 1;
                } else if(std::isdigit(yStr[0])){
                    Y = std::stoi(yStr);
                    z[i+1] = z[i] + Y;
                    B[n_i+newI-1][0] = Y;
                    B[n_i+newI-1][newI-1] = 1;
                }

            } else if (operation == "mul") {
                if (std::isdigit(xStr[0])) {
                    X = std::stoi(xStr);
                    z[i+1] = (z[i] * X) % p;
                    A[n_i+newI-1][newI-1] = X;
                    B[n_i+newI-1][newI-1] = 1;
                } else if (std::isdigit(yStr[0])) {
                    Y = std::stoi(yStr);
                    z[i+1] = (z[i] * Y) % p;
                    A[n_i+newI-1][newI-1] = 1;
                    B[n_i+newI-1][0] = Y;
                }
            }
        }
    }
    
    std::cout << "z[n]: ";
    for(int i = 0; i < n; i++){
        std::cout << z[i] << "\t";
    }
    std::cout << "\n";

    uint64_t H[n];
    uint64_t w, g_n;

    H[0] = 1;
    g_n = ((p - 1) / n) % p;
    w = power(2, g_n, p);
    for (int i = 1; i < n; i++)
    {
        H[i] = power(w, i, p);
    }
    std::cout << "H[n]: ";
    for(int i = 0; i < n; i++){
        std::cout << H[i] << "\t";
    }
    std::cout << "\n";

    uint64_t y, m, t, g_m;

    t = n_i + 1;
    m = (((power(n, 2, p) - n) / 2) - ((power(t, 2, p) - t) / 2)) % p;

    uint64_t K[m];
    K[0] = 1;
    g_m = ((p - 1) / m) % p;
    y = power(2, g_m, p);
    for (int i = 1; i < m; i++)
    {
        K[i] = power(y, i, p);
    }
    std::cout << "K[m]: ";
    for(int i = 0; i < m; i++){
        std::cout << K[i] << "\t";
    }
    std::cout << "\n";




    // std::vector<uint64_t> x(t), y(t);

    // Input the points
    std::vector<uint64_t> xL = {1, 43, 39};
    std::vector<uint64_t> yL = {42, 125, 135};

    // Calculate the interpolation
    auto [c2, c1, c0] = lagrange_interpolation(xL, yL, p);

    // Output the result
    std::cout << "row_A(x) = " << c2 << "x^2 + " << c1 << "x + " << c0 << " (mod " << p << ")" << std::endl;
    std::cout << "mod_inverse: " << mod_inverse(148, p) << "\n";

}

