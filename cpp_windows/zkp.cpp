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


// Function to multiply two polynomials
vector<uint64_t> multiplyPolynomials(const vector<uint64_t>& poly1, const vector<uint64_t>& poly2, uint64_t mod) {
    vector<uint64_t> result(poly1.size() + poly2.size() - 1, 0);

    for (size_t i = 0; i < poly1.size(); i++) {
        for (size_t j = 0; j < poly2.size(); j++) {
            result[i + j] = (result[i + j] + poly1[i] * poly2[j]) % mod;
        }
    }
    return result;
}

// Function to compute Lagrange basis polynomial L_i(x)
vector<uint64_t> LagrangePolynomial(int i, const vector<uint64_t>& x_values, uint64_t mod) {
    int n = x_values.size();
    vector<uint64_t> result = {1};  // Start with 1 for the polynomial (constant term)

    for (int j = 0; j < n; j++) {
        if (j != i) {
            vector<uint64_t> term = {static_cast<uint64_t>((mod - x_values[j]) % mod), 1};  // (x - x_j)
            uint64_t denominator = (x_values[i] + mod - x_values[j]) % mod;

            // Inverse of the denominator modulo mod
            uint64_t inv_denominator = 1;
            for (uint64_t d = 1; d < mod; d++) {
                if ((denominator * d) % mod == 1) {
                    inv_denominator = d;
                    break;
                }
            }

            // Multiply the result by (x - x_j) / (x_i - x_j)
            vector<uint64_t> temp = multiplyPolynomials(result, term, mod);
            for (uint64_t& coef : temp) {
                coef = (coef * inv_denominator) % mod;
            }
            result = temp;
        }
    }
    return result;
}

// Function to print polynomial
void PrintPolynomial(const vector<uint64_t>& poly, const string& name) {
    cout << name << "(x) = ";
    bool first = true;
    for (int i = poly.size() - 1; i >= 0; i--) {
        if (poly[i] != 0) {
            if (!first) {
                cout << " + ";
            }
            cout << poly[i] << "x^" << i;
            first = false;
        }
    }
    cout << endl;
}

void setupLagrangePolynomial (vector<uint64_t> x_values, vector<uint64_t> y_values, uint64_t mod) {
    // Automatically detect number of points
    int num_points = x_values.size();

    // Compute z-hatA(x) polynomial
    vector<uint64_t> z_hatA(1, 0);  // Start with a zero polynomial

    for (int i = 0; i < num_points; i++) {
        if (y_values[i] != 0) {  // Only process non-zero y-values
            vector<uint64_t> Li = LagrangePolynomial(i, x_values, mod);
            PrintPolynomial(Li, "L" + to_string(i + 1));

            // Multiply the L_i(x) by y_i and add to the final polynomial
            for (int j = 0; j < Li.size(); j++) {
                if (j >= z_hatA.size()) {
                    z_hatA.push_back(0);  // Ensure z_hatA is large enough to accommodate all terms
                }
                z_hatA[j] = (z_hatA[j] + y_values[i] * Li[j]) % mod;
            }
        }
    }
    // Print the final z-hatA(x) polynomial
    PrintPolynomial(z_hatA, "z-hatA");
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

    vector<vector<uint64_t>> Az(n, vector<uint64_t>(1, 0));
    vector<vector<uint64_t>> Bz(n, vector<uint64_t>(1, 0));
    vector<vector<uint64_t>> Cz(n, vector<uint64_t>(1, 0));
    // Matrix multiplication with modulo
    for (uint64_t i = 0; i < n; i++) {
        for (uint64_t j = 0; j < 1; j++) {
            for (uint64_t k = 0; k < n; k++) {
                Az[i][j] = (Az[i][j] + (A[i][k] * z[k]) % p) % p;
                Bz[i][j] = (Bz[i][j] + (B[i][k] * z[k]) % p) % p;
                Cz[i][j] = (Cz[i][j] + (C[i][k] * z[k]) % p) % p;
            }
        }
    }
    
    
    // Output the result
    cout << "Matrice Az under modulo " << p << " is:\n";
    for (uint64_t i = 0; i < n; i++) {
        cout << Az[i][0] << "\n";
    }
    cout << "Matrice Bz under modulo " << p << " is:" << "\n";
    for (uint64_t i = 0; i < n; i++) {
        cout << Bz[i][0] << "\n";
    }
    cout << "Matrice Cz under modulo " << p << " is:" << "\n";
    for (uint64_t i = 0; i < n; i++) {
        cout << Cz[i][0] << "\n";
    }

    int b = 2;
    
    
    // Given x-values and corresponding y-values
    vector<uint64_t> x_values;
    vector<uint64_t> y_values;

    
    uint64_t zA0 [n+b];
    uint64_t zA1 [n+b];
    cout << "zA(x) is:" << endl;
    for(uint64_t i = 0; i < n; i++) {
        zA0[i] = H[i];
        zA1[i] = Az[i][0];
        cout << "zA(" << zA0[i] << ")=" << zA1[i] << endl;
    }
    zA0[5] = 150;
    zA1[5] = 5;
    zA0[6] = 80;
    zA1[6] = 47;

    x_values.assign(zA0, zA0 + n+b);
    y_values.assign(zA1, zA1 + n+b);
    setupLagrangePolynomial(x_values, y_values, p);

    uint64_t zB0 [n+b];
    uint64_t zB1 [n+b];
    cout << "zB(x) is:" << endl;
    for(uint64_t i = 0; i < n; i++) {
        zB0[i] = H[i];
        zB1[i] = Bz[i][0];
        cout << "zB(" << zB0[i] << ")=" << zB1[i] << endl;
    }
    zB0[5] = zA0[5];
    zB1[5] = 15;
    zB0[6] = zA0[6];
    zB1[6] = 170;

    x_values.assign(zB0, zB0 + n+b);
    y_values.assign(zB1, zB1 + n+b);
    setupLagrangePolynomial(x_values, y_values, p);
    
    uint64_t zC0 [n+b];
    uint64_t zC1 [n+b];
    cout << "zC(x) is:" << endl;
    for(uint64_t i = 0; i < n; i++) {
        zC0[i] = H[i];
        zC1[i] = Cz[i][0];
        cout << "zC(" << zC0[i] << ")=" << zC1[i] << endl;
    }
    zC0[5] = zA0[5];
    zC1[5] = 1;
    zC0[6] = zA0[6];
    zC1[6] = 100;
    x_values.assign(zC0, zC0 + n+b);
    y_values.assign(zC1, zC1 + n+b);
    setupLagrangePolynomial(x_values, y_values, p);
    

    



}

