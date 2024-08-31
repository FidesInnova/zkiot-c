#include "zkp.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include <cstring>
#include <ctype.h> 
#include <cctype>

#include <tuple>
#include <unordered_map>
#include <string>
#include <map>

#include <random>
#include <set>
#include <algorithm>

#include <openssl/sha.h>
#include <iomanip>

#include "json.hpp"
// using json = nlohmann::json;
using ordered_json = nlohmann::ordered_json;
#include <regex>


using namespace std;


int64_t power(int64_t base, int64_t exponent, int64_t mod) {
    int64_t result = 1;
    base = base % mod;
    while (exponent > 0) {
        if (exponent % 2 == 1) {
            result = (result * base) % mod;
        }
        exponent = exponent >> 1;
        base = (base * base) % mod;
    }
    return result;
}

// Function to calculate the modular inverse using Extended Euclidean Algorithm
uint64_t mod_inverse(int64_t a, uint64_t mod) {
    int64_t t = 0, newt = 1;
    uint64_t r = mod, newr = a;

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


// Function to get the row indices of non-zero entries in matrix
vector<vector<int64_t>> getNonZeroRows(const vector<vector<int64_t>>& matrix) {
    uint64_t counter = 0;
    for (int64_t i = 0; i < matrix.size(); i++) {
        for (int64_t j = 0; j < matrix[i].size(); j++) {
            if (matrix[i][j] != 0) {
                counter++;
            }
        }
    }
    vector<vector<int64_t>> nonZeroRows(counter);
    counter = 0;
    for (int64_t i = 0; i < matrix.size(); i++) {
        for (int64_t j = 0; j < matrix[i].size(); j++) {
            if (matrix[i][j] != 0) {
                nonZeroRows[counter].push_back(i); // Storing row index
                nonZeroRows[counter].push_back(matrix[i][j]);
                counter++;
            }
        }
    }
    return nonZeroRows;
}


// Function to get the col indices of non-zero entries in matrix
vector<vector<int64_t>> getNonZeroCols(const vector<vector<int64_t>>& matrix) {
    uint64_t counter = 0;
    for (int64_t i = 0; i < matrix.size(); i++) {
        for (int64_t j = 0; j < matrix[i].size(); j++) {
            if (matrix[i][j] != 0) {
                counter++;
            }
        }
    }
    vector<vector<int64_t>> nonZeroCols(counter);
    counter = 0;
    for (int64_t i = 0; i < matrix.size(); i++) {
        for (int64_t j = 0; j < matrix[i].size(); j++) {
            if (matrix[i][j] != 0) {
                nonZeroCols[counter].push_back(j); // Storing col index
                nonZeroCols[counter].push_back(matrix[i][j]);
                counter++;
            }
        }
    }
    return nonZeroCols;
}


// Function to create the mapping
vector<vector<int64_t>> createMapping(const vector<int64_t>& K, const vector<int64_t>& H, const vector<vector<int64_t>>& nonZero) {
    vector<vector<int64_t>> row(K.size());
    for (uint64_t i = 0; i < nonZero.size(); i++) {
        row[i].push_back(K[i]);
        row[i].push_back(H[nonZero[i][0]]);
    }
    for (uint64_t i = nonZero.size(); i < row.size(); i++) {
        row[i].push_back(K[i]);
        row[i].push_back(H[i%H.size()]);
    }

    return row;
}


// Function to print the mapping
void printMapping(vector<vector<int64_t>>& row, const string& name) {
    for (uint64_t i = 0; i < row.size(); i++) {
        cout << name << "(" << row[i][0] << ") = " << row[i][1] << endl;
    }
    cout << endl;
}


// Function to create the val mapping
vector<vector<int64_t>> valMapping(const vector<int64_t>& K, const vector<int64_t>& H, vector<vector<int64_t>>& nonZeroRows, vector<vector<int64_t>>& nonZeroCols, uint64_t mod) {
    vector<vector<int64_t>> val(K.size());

    for (uint64_t i = 0; i < K.size(); i++) {
        if(i < nonZeroRows.size()) {
            val[i].push_back(K[i]);
            val[i].push_back((nonZeroRows[i][1]*mod_inverse(((H.size()*power(H[nonZeroRows[i][0]], H.size()-1, mod))%mod) * ((H.size()*power(H[nonZeroCols[i][0]], H.size()-1, mod))%mod), mod))%mod);
        }
        else {
            val[i].push_back(K[i]);
            val[i].push_back(0);
        }
    }

    return val;
}


// Function to expand polynomials given the roots
std::vector<int64_t> expandPolynomials(const std::vector<int64_t>& roots, int64_t mod) {
    std::vector<int64_t> result = {1}; // Start with the polynomial "1"

    for (int64_t root : roots) {
        // Multiply the current result polynomial by (x - root)
        std::vector<int64_t> temp(result.size() + 1, 0);
        for (size_t i = 0; i < result.size(); i++) {
            temp[i] += result[i]; // x^n term
            temp[i] %= mod;
            temp[i + 1] -= result[i] * root; // -root * x^(n-1) term
            temp[i + 1] %= mod;
        }
        result = temp;
    }
    std::reverse(result.begin(), result.end());
    return result;
}

// Function to parse a polynomial string and return it as a vector
vector<int64_t> parsePolynomial(const string& poly) {
    map<int64_t, int64_t> terms; // Use int64_t to match your usage
    istringstream iss(poly);
    string term;

    // Improved parsing logic
    while (getline(iss, term, '+')) {
        int64_t coef = 0;
        int64_t exp = 0;
        size_t x_pos = term.find("x^");

        if (x_pos != string::npos) {
            coef = stoll(term.substr(0, x_pos));  // Get the coefficient before "x^"
            exp = stoi(term.substr(x_pos + 2));   // Get the exponent after "x^"
        } else {
            // If there's no "x^", it's a constant term
            coef = stoll(term);
            exp = 0;
        }
        
        terms[exp] = coef;
    }

    // Determine the degree of the polynomial
    int64_t degree = terms.rbegin()->first;

    // Prepare the vector to represent the polynomial
    vector<int64_t> coefficients(degree + 1, 0);
    for (const auto& term : terms) {
        coefficients[term.first] = term.second;
    }

    return coefficients;
}


// Function to add two polynomials
std::vector<int64_t> addPolynomials(const std::vector<int64_t>& poly1, const std::vector<int64_t>& poly2, int64_t mod) {
    // Determine the size of the result polynomial (the max size of the two input polynomials)
    size_t maxSize = std::max(poly1.size(), poly2.size());
    std::vector<int64_t> result(maxSize, 0);

    // Add the polynomials
    for (size_t i = 0; i < maxSize; ++i) {
        if (i < poly1.size()) {
            result[i] = (result[i] + poly1[i]) % mod;
        }
        if (i < poly2.size()) {
            result[i] = (result[i] + poly2[i]) % mod;
        }
        // Ensure the result is non-negative
        if (result[i] < 0) {
            result[i] += mod;
        }
    }

    return result;
}


// Function to subtract one polynomial from another
vector<int64_t> subtractPolynomials(const vector<int64_t>& poly1, const vector<int64_t>& poly2, int64_t mod) {
    size_t maxSize = std::max(poly1.size(), poly2.size());
    vector<int64_t> buf1(maxSize, 0);
    vector<int64_t> buf2(maxSize, 0);
    vector<int64_t> result(maxSize, 0);
    buf1 = poly1;
    buf2 = poly2;
    for (size_t i = 0; i < maxSize; i++) {
        result[i] = (buf1[i] - buf2[i]) % mod;
        if (result[i] < 0) {
            result[i] += mod;  // Handle negative results
        }
    }
    return result;
}


// Function to multiply two polynomials
vector<int64_t> multiplyPolynomials(const vector<int64_t>& poly1, const vector<int64_t>& poly2, int64_t mod) {
    vector<int64_t> result(poly1.size() + poly2.size() - 1, 0);

    for (size_t i = 0; i < poly1.size(); i++) {
        for (size_t j = 0; j < poly2.size(); j++) {
            result[i + j] = (result[i + j] + poly1[i] * poly2[j]) % mod;
            if (result[i + j] < 0) {
                result[i + j] += mod;
            }
        }
    }
    return result;
}




int64_t generateRandomNumber(const std::vector<int64_t>& H, int64_t mod) {
    std::mt19937_64 rng(std::random_device{}());  // Use random_device to seed the generator
    std::uniform_int_distribution<int64_t> dist(0, mod - 1);
    
    int64_t randomNumber;
    do {
        randomNumber = dist(rng);
    } while (std::find(H.begin(), H.end(), randomNumber) != H.end());

    return randomNumber;
}

// Function to convert int64_t to a byte array
void int64ToByteArray(int64_t value, unsigned char* byteArray) {
    for (size_t i = 0; i < sizeof(value); ++i) {
        byteArray[i] = static_cast<unsigned char>(value >> (i * 8));
    }
}


// Function to divide a polynomial by another polynomial
vector<vector<int64_t>> dividePolynomials(const vector<int64_t>& dividend, const vector<int64_t>& divisor, int64_t mod) {
    vector<int64_t> quotient;
    vector<int64_t> remainder = dividend;
    vector<vector<int64_t>> result;

    int64_t n = remainder.size();
    int64_t m = divisor.size();


    // If the divisor is larger than the dividend, the quotient is zero
    if (m > n) {
        quotient.push_back(0);
        result.push_back(quotient);
        result.push_back(remainder);
        return result;
    }

    // Main loop for polynomial division
    while (remainder.size() >= divisor.size()) {
        // Degree of the current term in the quotient
        int64_t degreeDiff = remainder.size() - divisor.size();
        int64_t coef = remainder.back() * mod_inverse(divisor.back(), mod) % mod;
        vector<int64_t> term(degreeDiff + 1, 0);
        term[degreeDiff] = coef;

        quotient.insert(quotient.begin(), coef);
        vector<int64_t> termTimesDivisor = multiplyPolynomials(term, divisor, mod);
        termTimesDivisor.resize(remainder.size(), 0);
        remainder = subtractPolynomials(remainder, termTimesDivisor, mod);

        // Remove leading zeros
        while (remainder.size() > 0 && remainder.back() == 0) {
            remainder.pop_back();
        }
    }
    

    if (quotient.empty()) {
        quotient.push_back(0);
    }
    
    result.push_back(quotient);
    result.push_back(remainder);
    return result;
}




// Function to calculate Polynomial r(α,x) = (alpha^n - x^n) / (alpha - x)
std::vector<int64_t> calculatePolynomial_r_alpha_x(int64_t alpha, int64_t n, int64_t mod) {
    std::vector<int64_t> P(n, 0);

    // Calculate each term of the polynomial P(x)
    int64_t currentPowerOfAlpha = 1; // alpha^0
    for (int64_t i = 0; i < n; ++i) {
        P[n - 1 - i] = currentPowerOfAlpha; // alpha^(n-1-i)
        currentPowerOfAlpha = (currentPowerOfAlpha * alpha) % mod;
    }

    return P;
}


// Function to calculate Polynomial r(α,x) = (alpha^n - x^n) / (alpha - x)
int64_t calculatePolynomial_r_alpha_k(int64_t alpha, int64_t k, int64_t n, int64_t mod) {
    int64_t result = 1;
    result = (power(alpha, n, mod)-power(k, n, mod)) % mod;
    int64_t buff = (alpha - k) % mod;
    if (buff < 0) {
        buff += mod;
    }
    result *= mod_inverse(buff, mod);;
    result %= mod;
    return result;
}


// Function to generate a random polynomial
vector<int64_t> generateRandomPolynomial(size_t numTerms, size_t maxDegree, int64_t mod) {
    if (numTerms > maxDegree + 1) {
        throw invalid_argument("Number of terms cannot be greater than the number of possible degrees.");
    }

    vector<int64_t> polynomial(maxDegree + 1, 0); // Initialize polynomial with zeros

    // Generate random indices for the non-zero terms
    set<size_t> indices;
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<size_t> dis(0, maxDegree);

    while (indices.size() < numTerms) {
        indices.insert(dis(gen));
    }

    // Fill the polynomial with random values at the chosen indices
    for (size_t index : indices) {
        polynomial[index] = dis(gen) % mod;
        if (polynomial[index] < 0) polynomial[index] += mod;
    }

    return polynomial;
}


// Function to create a polynomial for (x - root)
vector<int64_t> createLinearPolynomial(int64_t root) {
    return {-root, 1};  // Represents (x - root)
}


// Function to compute Lagrange basis polynomial L_i(x)
vector<int64_t> LagrangePolynomial(int64_t i, const std::vector<int64_t>& x_values, int64_t mod) {
    int64_t n = x_values.size();
    vector<int64_t> result = {1};  // Start with 1 for the polynomial (constant term)

    for (int64_t j = 0; j < n; j++) {
        if (j != i) {
            vector<int64_t> term = {static_cast<int64_t>((mod - x_values[j]) % mod), 1};  // (x - x_j)
            int64_t denominator = (x_values[i] + mod - x_values[j]) % mod;

            // Inverse of the denominator modulo mod
            int64_t inv_denominator = 1;
            for (int64_t d = 1; d < mod; d++) {
                if ((denominator * d) % mod == 1) {
                    inv_denominator = d;
                    break;
                }
            }

            // Multiply the result by (x - x_j) / (x_i - x_j)
            vector<int64_t> temp = multiplyPolynomials(result, term, mod);
            for (int64_t& coef : temp) {
                coef = (coef * inv_denominator) % mod;
            }
            result = temp;
        }
    }
    return result;
}


// Function to print polynomial
void PrintPolynomial(const vector<int64_t>& poly, const string& name) {
    cout << name;
    bool first = true;
    for (int64_t i = poly.size() -1 ; i >= 0; i--) {
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

// Global map to hold the dynamically created global variables
std::unordered_map<const char*, std::string> LagrangeOutput;
// Function to store polynomial
void storePolynomial(const std::vector<int64_t>& poly, const char* name) {
    bool first = true;
    std::string buffer; 
    for (int64_t i = poly.size() - 1; i >= 0; i--) {
        if (poly[i] != 0) {
            if (!first) {
                buffer += "+";
            }
            buffer += to_string(poly[i]) + "x^" + to_string(i);
            first = false;
        }
    }
    LagrangeOutput[name] = buffer;
    cout << name << " = " << LagrangeOutput[name] << endl;
}


void setupLagrangePolynomial (const std::vector<int64_t> x_values, const std::vector<int64_t> y_values, int64_t mod, const char* name) {
    // Automatically detect number of points
    int64_t num_points = x_values.size();
    
    // Compute polynomial(x) polynomial
    vector<int64_t> polynomial(1, 0);  // Start with a zero polynomial

    for (int64_t i = 0; i < num_points; i++) {
        if (y_values[i] != 0) {  // Only process non-zero y-values
            vector<int64_t> Li = LagrangePolynomial(i, x_values, mod);
            PrintPolynomial(Li, "L" + to_string(i + 1) + ": ");

            // Multiply the L_i(x) by y_i and add to the final polynomial
            for (int64_t j = 0; j < Li.size(); j++) {
                if (j >= polynomial.size()) {
                    polynomial.push_back(0);  // Ensure polynomial is large enough to accommodate all terms
                }
                polynomial[j] = (polynomial[j] + y_values[i] * Li[j]) % mod;
            }
        }
    }
    // Print the final polynomial(x) polynomial
    storePolynomial(polynomial, name);
}


void ZKP::setup(int64_t g, int64_t d, int64_t l, int64_t p, std::vector<int64_t>& pp) {
    pp.clear();
    int64_t pMinusOne = p - 1;
    int64_t current_exponent = 1;
    int64_t newPower = d % pMinusOne;
    int64_t res = 0;

    res = power(g, newPower, p);
    pp.push_back(res);
    for (int64_t i = 1; i < l; ++i) {
        res = power(res, newPower, p);
        pp.push_back(res);
    }
    //pp=(2,66,83,91,96,24,2,66,83)
}


// Function to parse the polynomial string and evaluate it
int64_t evaluatePolynomial(const vector<int64_t>& polynomial, int64_t x, int64_t mod) {
    int64_t result = 0;

    for (size_t i = 0; i < polynomial.size(); i++) {
        
        int64_t termValue = polynomial[i] * power(x, i, mod) % mod;
        
        if (termValue < 0) {
            termValue += mod;
        }
        result += termValue;
        result %= mod;
        if (result < 0) {
            result += mod;
        }
    }
    return result;
}



// Function to compute the sum of polynomial evaluations at multiple points
int64_t sumOfEvaluations(const vector<int64_t>& poly, const vector<int64_t>& points, int64_t mod) {
    int64_t totalSum = 0;

    for (int64_t point : points) {
        int64_t evlp = evaluatePolynomial(poly, point, mod);
        totalSum = (totalSum + evlp) % mod;
    }
    totalSum %= mod;

    if (totalSum < 0) totalSum += mod;

    return totalSum;
}


vector<int64_t> multiplyPolynomialByNumber(const vector<int64_t>& H, int64_t h, int64_t mod) {
    vector<int64_t> result(H.size(), 1);  // Use long long to avoid overflow during multiplication

    for (int64_t i = 0; i < H.size(); i++) {
        result [i] = (H[i] * h) % mod;

        // // If result becomes negative, convert it to a positive equivalent under modulo
        if (result[i] < 0) {
            result[i] += mod;
        }
    }
    return result;
}

// Helper function to initialize a zero matrix of given size
std::vector<std::vector<int64_t>> initializeZeroMatrix(int64_t size) {
    return std::vector<std::vector<int64_t>>(size, std::vector<int64_t>(size, 0ll));
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


void ZKP::createMat(const std::string& filename, std::vector<std::vector<int64_t>>& A, std::vector<std::vector<int64_t>>& B, std::vector<std::vector<int64_t>>& C, int64_t p) {
    // Reading the instructions from the file
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("Error opening file");
    }

    std::string line;
    std::vector<std::string> instructions;
    std::unordered_map<int64_t, int64_t> inputs; // Store which registers are inputs
    std::cout << "Reading " << filename << "\n";
    while (std::getline(file, line)) {
        // std::cout << line << "\n";
        instructions.push_back(line);
    }
    // std::cout << "End of the file\n";
    file.close();

    // Number of gates and inputs
    int64_t n_g = 0;
    int64_t n_i = 0;
    
    // Parse instructions to determine n_g and n_i
    for (const auto& instr : instructions) {
        std::stringstream ss(instr);
        std::string operation, rStr, xStr, yStr;
        
        ss >> operation >> rStr;

        if (operation == "li") {
            ss >> xStr;
            n_i++;
            int64_t R = std::stoi(rStr.substr(1));
            int64_t X = std::stoi(xStr);
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
    int64_t n = n_g + n_i + 1;
    std::cout << "Matrix order: " << n << "\n";

    // Initialize matrices A, B, C
    A = initializeZeroMatrix(n);
    B = initializeZeroMatrix(n);
    C = initializeZeroMatrix(n);


    // Fill matrices based on the instructions
    int64_t gateIndex = 0;
    int64_t z[n+1];
    z[0] = 1;

    for (int64_t i = 0; i < n-1; i++) {
        std::stringstream ss(instructions[i]);
        std::string operation, rStr, xStr, yStr;
        ss >> operation >> rStr;
        int64_t R = std::stoi(rStr.substr(1));

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

        int64_t X,Y;
        if (operation == "li") {
            X = std::stoi(xStr);
            z[i+1] = X % p;
        } else {
            ss >> xStr >> yStr;
            
            gateIndex++;
            int64_t newI = gateIndex + n_i;
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
    
    cout << "z[n]: ";
    for (int64_t i = 0; i < n; i++){
        cout << z[i] << " ";
    }
    cout << endl;

    vector <int64_t> H;
    int64_t w, g_n;

    H.push_back(1);
    // H[0] = 1;
    g_n = ((p - 1) / n) % p;
    w = power(2, g_n, p);
    for (int64_t i = 1; i < n; i++)
    {
        // H[i] = power(w, i, p);
        H.push_back(power(w, i, p));
    }
    cout << "H[n]: ";
    for (int64_t i = 0; i < n; i++){
        cout << H[i] << " ";
    }
    cout << endl;

    int64_t y, m, t, g_m;

    t = n_i + 1;
    m = (((power(n, 2, p) - n) / 2) - ((power(t, 2, p) - t) / 2)) % p;

    vector<int64_t> K;
    K.push_back(1);
    g_m = ((p - 1) / m) % p;
    y = power(2, g_m, p);
    for (int64_t i = 1; i < m; i++)
    {
        K.push_back(power(y, i, p));
    }
    cout << "K[m]: ";
    for (int64_t i = 0; i < m; i++){
        cout << K[i] << " ";
    }
    cout << endl;

    vector<vector<int64_t>> Az(n, vector<int64_t>(1, 0));
    vector<vector<int64_t>> Bz(n, vector<int64_t>(1, 0));
    vector<vector<int64_t>> Cz(n, vector<int64_t>(1, 0));
    // Matrix multiplication with modulo
    for (int64_t i = 0; i < n; i++) {
        for (int64_t j = 0; j < 1; j++) {
            for (int64_t k = 0; k < n; k++) {
                Az[i][j] = (Az[i][j] + (A[i][k] * z[k]) % p) % p;
                Bz[i][j] = (Bz[i][j] + (B[i][k] * z[k]) % p) % p;
                Cz[i][j] = (Cz[i][j] + (C[i][k] * z[k]) % p) % p;
            }
        }
    }
    
    
    // Output the result
    cout << "Matrice Az under modulo " << p << " is: ";
    for (int64_t i = 0; i < n; i++) {
        cout << Az[i][0] << " ";
    }
    cout << endl;
    cout << "Matrice Bz under modulo " << p << " is:" << " ";
    for (int64_t i = 0; i < n; i++) {
        cout << Bz[i][0] << " ";
    }
    cout << endl;
    cout << "Matrice Cz under modulo " << p << " is:" << " ";
    for (int64_t i = 0; i < n; i++) {
        cout << Cz[i][0] << " ";
    }
    cout << endl;

    int64_t b = 2;
    
    
    // Given x-values and corresponding y-values
    vector<int64_t> x_values;
    vector<int64_t> y_values;

    
    int64_t zA0 [n+b];
    int64_t zA1 [n+b];
    cout << "zA(x) is:" << endl;
    for (int64_t i = 0; i < n; i++) {
        zA0[i] = H[i];
        zA1[i] = Az[i][0];
        cout << "zA(" << zA0[i] << ") = " << zA1[i] << endl;
    }
    zA0[5] = 150;
    zA1[5] = 5;
    zA0[6] = 80;
    zA1[6] = 47;
    x_values.assign(zA0, zA0 + n+b);
    y_values.assign(zA1, zA1 + n+b);
    setupLagrangePolynomial(x_values, y_values, p, "z_hatA(x)");

    int64_t zB0 [n+b];
    int64_t zB1 [n+b];
    cout << "zB(x) is:" << endl;
    for (int64_t i = 0; i < n; i++) {
        zB0[i] = H[i];
        zB1[i] = Bz[i][0];
        cout << "zB(" << zB0[i] << ") = " << zB1[i] << endl;
    }
    zB0[5] = zA0[5];
    zB1[5] = 15;
    zB0[6] = zA0[6];
    zB1[6] = 170;

    x_values.assign(zB0, zB0 + n+b);
    y_values.assign(zB1, zB1 + n+b);
    setupLagrangePolynomial(x_values, y_values, p, "z_hatB(x)");
    
    int64_t zC0 [n+b];
    int64_t zC1 [n+b];
    cout << "zC(x) is:" << endl;
    for (int64_t i = 0; i < n; i++) {
        zC0[i] = H[i];
        zC1[i] = Cz[i][0];
        cout << "zC(" << zC0[i] << ") = " << zC1[i] << endl;
    }
    zC0[5] = zA0[5];
    zC1[5] = 1;
    zC0[6] = zA0[6];
    zC1[6] = 100;
    x_values.assign(zC0, zC0 + n+b);
    y_values.assign(zC1, zC1 + n+b);
    setupLagrangePolynomial(x_values, y_values, p, "z_hatC(x)");
    
    
    int64_t zero_to_t_for_H[t];
    int64_t t_to_n_for_H[n-t];
    int64_t zero_to_t_for_z[t];
    int64_t t_to_n_for_z[n-t];

    for (int64_t i = 0; i < t; i++) {
        zero_to_t_for_H[i] = H[i];
        zero_to_t_for_z[i] = z[i];
    }
    for (int64_t i = 0; i < n-t; i++) {
        t_to_n_for_H[i] = H[i+t];
        t_to_n_for_z[i] = z[i+t];
    }

    x_values.assign(zero_to_t_for_H, zero_to_t_for_H + t);
    y_values.assign(zero_to_t_for_z, zero_to_t_for_z + t);

    setupLagrangePolynomial(x_values, y_values, p, "x_hat(h)");


    cout << "w_bar(h) is:" << endl;
    vector<int64_t> w_bar(n-t + b);
    x_values.assign(zero_to_t_for_H, zero_to_t_for_H + t);
    vector<int64_t> polyX_HAT_H = parsePolynomial(LagrangeOutput["x_hat(h)"]);
    vector<int64_t> w_bar_numerator (n-t, 1);
    vector<int64_t> w_bar_denominator (n-t, 1);
    for (int64_t i = 0; i < n-t; i++){
        w_bar_numerator[i] = (t_to_n_for_z[i] - (evaluatePolynomial(polyX_HAT_H, t_to_n_for_H[i], p)));
        
        for (int64_t j = 0; j < x_values.size(); j++) {
            w_bar_denominator[i] *= (t_to_n_for_H[i] - x_values[j]);
            // Apply modulus to keep the number within the bounds
            w_bar_denominator[i] %= p;

            // // If result becomes negative, convert it to a positive equivalent under modulo
            if (w_bar_denominator[i] < 0) {
                w_bar_denominator[i] += p;
            }
        }
        w_bar_denominator[i] = mod_inverse(w_bar_denominator[i], p);
        w_bar[i] = (w_bar_numerator[i] * w_bar_denominator[i]) %p;
        
        if (w_bar[i] < 0) {
            w_bar[i] += p;
        }
        cout << "w_bar(" << t_to_n_for_H[i] << ") = " << w_bar[i] << endl;
    }

    int64_t w_hat0[n-t+b];
    int64_t w_hat1[n-t+b];

    for (int64_t i = 0; i < n-t; i++) {
        w_hat0[i] = t_to_n_for_H[i];
        w_hat1[i] = w_bar[i];
    }

    w_hat0[3] = 150;
    w_hat1[3] = 42;
    w_hat0[4] = 80;
    w_hat1[4] = 180;

    x_values.assign(w_hat0, w_hat0 + n-t+b);
    y_values.assign(w_hat1, w_hat1 + n-t+b);
    setupLagrangePolynomial(x_values, y_values, p, "w_hat(x)");

    // Parse the polynomials
    vector<int64_t> z_hatA = parsePolynomial(LagrangeOutput["z_hatA(x)"]);
    vector<int64_t> z_hatB = parsePolynomial(LagrangeOutput["z_hatB(x)"]);
    vector<int64_t> z_hatC = parsePolynomial(LagrangeOutput["z_hatC(x)"]);

    // Multiply A and B
    vector<int64_t> productAB = multiplyPolynomials(z_hatA, z_hatB, p);
    // Subtract C from the product of A and B
    vector<int64_t> zAzB_zC = subtractPolynomials(productAB, z_hatC, p);

    storePolynomial(zAzB_zC, "zA(x)zB(x)-zC(x)");



    // Start with the first polynomial (x - H)
    vector<int64_t> vH_x = createLinearPolynomial(H[0]);

    // Multiply (x - H) for all other H
    for (size_t i = 1; i < n; i++) {
        vector<int64_t> nextPoly = createLinearPolynomial(H[i]);
        vH_x = multiplyPolynomials(vH_x, nextPoly, p);
    }
    storePolynomial(vH_x, "vH(x)");


    // Example division: Dividing the product of zAzB_zC by vH_x
    vector<int64_t> h_0_x = dividePolynomials(zAzB_zC, vH_x, p)[0];

    storePolynomial(h_0_x, "h0(x)");

    // vector<int64_t> s_x = generateRandomPolynomial(n, (2*n)+b-1, p);
    vector<int64_t> s_x = {115, 3, 0, 0, 20, 1, 0, 17, 101, 0, 5};
    storePolynomial(s_x, "s(x)");


    // x_values.assign(H, H + n);
    int64_t sigma1 = sumOfEvaluations(s_x, H, p);
    cout << "sigma1(H) = " << sigma1 << endl;


    // Create a random device and a Mersenne Twister random number generator
    std::random_device rd;
    std::mt19937_64 gen(rd());
    // Define a uniform distribution for uint64_t numbers in the range [0, mod-1]
    std::uniform_int_distribution<uint64_t> dis(0, p - 1);

    // Generate random numbers with the modulus applied
    // int64_t alpha = dis(gen);
    // int64_t etaA = dis(gen);
    // int64_t etaB = dis(gen);
    // int64_t etaC = dis(gen);

    int64_t alpha = 10;
    int64_t etaA = 2;
    int64_t etaB = 30;
    int64_t etaC = 100;

    // Print the generated numbers
    cout << "alpha: " << alpha << endl;
    cout << "etaA: " << etaA << endl;
    cout << "etaB: " << etaB << endl;
    cout << "etaC: " << etaC << endl;

    vector<int64_t> etaA_z_hatA_x(n+b);
    vector<int64_t> etaB_z_hatB_x(n+b);
    vector<int64_t> etaC_z_hatC_x(n+b);
    
    etaA_z_hatA_x = multiplyPolynomialByNumber(z_hatA, etaA, p);
    etaB_z_hatB_x = multiplyPolynomialByNumber(z_hatB, etaB, p);
    etaC_z_hatC_x = multiplyPolynomialByNumber(z_hatC, etaC, p);

    vector<int64_t> Sum_M_eta_M_z_hat_M_x = addPolynomials(addPolynomials(etaA_z_hatA_x, etaB_z_hatB_x, p), etaC_z_hatC_x, p);
    storePolynomial(Sum_M_eta_M_z_hat_M_x, "Sum_M_z_hatM(x)");

    vector<int64_t> r_alpha_x = calculatePolynomial_r_alpha_x(alpha, n, p);
    storePolynomial(r_alpha_x, "r(alpha, x)");

    vector<int64_t> r_Sum_x = multiplyPolynomials(r_alpha_x, Sum_M_eta_M_z_hat_M_x, p);
    storePolynomial(r_Sum_x, "r(alpha, x)Sum_M_z_hatM(x)");

    vector<int64_t> w_hat_x = parsePolynomial(LagrangeOutput["w_hat(x)"]);
    
    
    vector<int64_t> v_H;
    v_H.assign(zero_to_t_for_H, zero_to_t_for_H + t);
    v_H = expandPolynomials(v_H, p);
    storePolynomial(v_H, "v_H");
    vector<int64_t> z_hat_x = addPolynomials(multiplyPolynomials(w_hat_x, v_H, p), polyX_HAT_H, p);
    storePolynomial(z_hat_x, "z_hat(x)");



    // Get the non-zero rows from matrix A
    vector<vector<int64_t>> nonZeroRowsA = getNonZeroRows(A);
    // Create the rowA mapping
    vector<vector<int64_t>> rowA;
    rowA = createMapping(K, H, nonZeroRowsA);
    // Print the rowA mapping
    rowA[3][1] = 1;
    rowA[4][1] = 135;
    rowA[5][1] = 125;
    rowA[6][1] = 59;
    rowA[7][1] = 42;
    rowA[8][1] = 1;
    printMapping(rowA, "row_A");

    // Get the non-zero cols from matrix A
    vector<vector<int64_t>> nonZeroColsA = getNonZeroCols(A);
    // Create the colA mapping
    vector<vector<int64_t>> colA;
    colA = createMapping(K, H, nonZeroColsA);
    // Print the colA mapping
    colA[3][1] = 42;
    colA[4][1] = 1;
    colA[5][1] = 135;
    colA[6][1] = 125;
    colA[7][1] = 59;
    colA[8][1] = 42;
    printMapping(colA, "col_A");
    
    vector<vector<int64_t>> valA = valMapping(K, H, nonZeroRowsA, nonZeroColsA, p);
    printMapping(valA, "val_A");

    vector<int64_t> A_hat(H.size(), 0);
    for (uint64_t i = 0; i < nonZeroRowsA.size(); i++) {
        vector<int64_t> buff_n = calculatePolynomial_r_alpha_x(rowA[i][1], H.size(), p);
        int64_t eval = evaluatePolynomial(buff_n, rowA[i][1], p);
        vector<int64_t> buff = calculatePolynomial_r_alpha_x(colA[i][1], H.size(), p);
        eval = (eval *valA[i][1]) % p;
        if (eval < 0) {
            eval += p;
        }
        buff = multiplyPolynomialByNumber(buff, eval, p);
        buff = multiplyPolynomialByNumber(buff, calculatePolynomial_r_alpha_k(alpha, rowA[i][1], H.size(), p), p);
        PrintPolynomial(buff, "A_hat(x): ");
        if(i > 0) {
            A_hat = addPolynomials(A_hat, buff, p);
        }
        else {
            A_hat = buff;
        }
    }
    storePolynomial(A_hat, "A_hat");



    vector<vector<int64_t>> nonZeroRowsB = getNonZeroRows(B);
    vector<vector<int64_t>> rowB;
    rowB = createMapping(K, H, nonZeroRowsB);
    rowB[4][1] = 59;
    rowB[5][1] = 1;
    rowB[6][1] = 42;
    rowB[7][1] = 135;
    rowB[8][1] = 59;
    printMapping(rowB, "row_B");

    vector<vector<int64_t>> nonZeroColsB = getNonZeroCols(B);
    vector<vector<int64_t>> colB;
    colB = createMapping(K, H, nonZeroColsB);
    colB[4][1] = 59;
    colB[5][1] = 42;
    colB[6][1] = 125;
    colB[7][1] = 1;
    colB[8][1] = 135;
    printMapping(colB, "col_B");

    vector<vector<int64_t>> valB = valMapping(K, H, nonZeroRowsB, nonZeroColsB, p);
    printMapping(valB, "val_B");

    vector<int64_t> B_hat(H.size(), 0);
    for (uint64_t i = 0; i < nonZeroRowsB.size(); i++) {
        vector<int64_t> buff_n = calculatePolynomial_r_alpha_x(rowB[i][1], H.size(), p);
        int64_t eval = evaluatePolynomial(buff_n, rowB[i][1], p);
        vector<int64_t> buff = calculatePolynomial_r_alpha_x(colB[i][1], H.size(), p);
        eval = (eval *valB[i][1]) % p;
        if (eval < 0) {
            eval += p;
        }
        buff = multiplyPolynomialByNumber(buff, eval, p);
        buff = multiplyPolynomialByNumber(buff, calculatePolynomial_r_alpha_k(alpha, rowB[i][1], H.size(), p), p);
        PrintPolynomial(buff, "B_hat(x): ");
        if(i > 0) {
            B_hat = addPolynomials(B_hat, buff, p);
        }
        else {
            B_hat = buff;
        }
    }
    storePolynomial(B_hat, "B_hat");




    vector<vector<int64_t>> nonZeroRowsC = getNonZeroRows(C);
    vector<vector<int64_t>> rowC;
    rowC = createMapping(K, H, nonZeroRowsC);
    rowC[3][1] = 1;
    rowC[4][1] = 59;
    rowC[5][1] = 125;
    rowC[6][1] = 1;
    rowC[7][1] = 135;
    rowC[8][1] = 42;
    printMapping(rowC, "row_C");

    vector<vector<int64_t>> nonZeroColsC = getNonZeroCols(C);
    vector<vector<int64_t>> colC;
    colC = createMapping(K, H, nonZeroColsC);
    colC[3][1] = 125;
    colC[4][1] = 59;
    colC[5][1] = 1;
    colC[6][1] = 1;
    colC[7][1] = 42;
    colC[8][1] = 59;
    printMapping(colC, "col_C");
    
    vector<vector<int64_t>> valC = valMapping(K, H, nonZeroRowsC, nonZeroColsC, p);
    printMapping(valC, "val_C");

    vector<int64_t> C_hat(H.size(), 0);
    for (uint64_t i = 0; i < nonZeroRowsC.size(); i++) {
        vector<int64_t> buff_n = calculatePolynomial_r_alpha_x(rowC[i][1], H.size(), p);
        int64_t eval = evaluatePolynomial(buff_n, rowC[i][1], p);
        vector<int64_t> buff = calculatePolynomial_r_alpha_x(colC[i][1], H.size(), p);
        eval = (eval *valC[i][1]) % p;
        if (eval < 0) {
            eval += p;
        }
        buff = multiplyPolynomialByNumber(buff, eval, p);
        buff = multiplyPolynomialByNumber(buff, calculatePolynomial_r_alpha_k(alpha, rowC[i][1], H.size(), p), p);
        PrintPolynomial(buff, "C_hat(x): ");
        if(i > 0) {
            C_hat = addPolynomials(C_hat, buff, p);
        }
        else {
            C_hat = buff;
        }
    }
    storePolynomial(C_hat, "C_hat");



    vector<int64_t> eta_A_hat = multiplyPolynomialByNumber(A_hat, etaA, p);
    vector<int64_t> eta_B_hat = multiplyPolynomialByNumber(B_hat, etaB, p);
    vector<int64_t> eta_C_hat = multiplyPolynomialByNumber(C_hat, etaC, p);
    
    PrintPolynomial(eta_A_hat, "eta_A_hat: ");
    PrintPolynomial(eta_B_hat, "eta_B_hat: ");
    PrintPolynomial(eta_C_hat, "eta_C_hat: ");

    vector<int64_t> Sum_M_eta_M_r_M_alpha_x = addPolynomials(addPolynomials(eta_A_hat, eta_B_hat, p), eta_C_hat, p);

    storePolynomial(Sum_M_eta_M_r_M_alpha_x, "Sum_M_eta_M_r_M(alpha ,x)");


    vector<int64_t> Sum_M_eta_M_r_M_alpha_x_z_hat_x = multiplyPolynomials(Sum_M_eta_M_r_M_alpha_x, z_hat_x, p);


    storePolynomial(Sum_M_eta_M_r_M_alpha_x_z_hat_x, "Sum_M_eta_M_r_M(alpha ,x)z-hat(x)");

    vector<int64_t> Sum_check_protocol = addPolynomials(s_x, (subtractPolynomials(r_Sum_x, Sum_M_eta_M_r_M_alpha_x_z_hat_x, p)), p);
    storePolynomial(Sum_check_protocol, "Sum_check_protocol");

    vector<int64_t> h_1_x = dividePolynomials(Sum_check_protocol, vH_x, p)[0];
    storePolynomial(h_1_x, "h1(x)");

    vector<int64_t> g_1_x = dividePolynomials(Sum_check_protocol, vH_x, p)[1];
    g_1_x.erase(g_1_x.begin());
    storePolynomial(g_1_x, "g1(x)");

    // int64_t beta1 = generateRandomNumber(H, p);
    int64_t beta1 = 22;
    cout << "beta1 = " << beta1 << endl;


    int64_t sigma2 = 0;
    sigma2 = ((etaA * evaluatePolynomial(A_hat, beta1, p)) % p + (etaB * evaluatePolynomial(B_hat, beta1, p)) % p + (etaC * evaluatePolynomial(C_hat, beta1, p)) % p) % p;
    cout << "sigma2 = " << sigma2 << endl;



    vector<int64_t> A_hat_M_hat(H.size(), 0);
    vector<int64_t> B_hat_M_hat(H.size(), 0);
    vector<int64_t> C_hat_M_hat(H.size(), 0);
    for (uint64_t i = 0; i < nonZeroRowsA.size(); i++) {
        vector<int64_t> buff_nA = calculatePolynomial_r_alpha_x(colA[i][1], H.size(), p);
        int64_t evalA = evaluatePolynomial(buff_nA, beta1, p);
        vector<int64_t> buffA = calculatePolynomial_r_alpha_x(rowA[i][1], H.size(), p);
        evalA = (evalA *valA[i][1]) % p;
        if (evalA < 0) {
            evalA += p;
        }
        buffA = multiplyPolynomialByNumber(buffA, evalA, p);

        if(i > 0) {
            A_hat_M_hat = addPolynomials(A_hat_M_hat, buffA, p);
        }
        else {
            A_hat_M_hat = buffA;
        }
    }
    for (uint64_t i = 0; i < nonZeroRowsB.size(); i++) {
        vector<int64_t> buff_nB = calculatePolynomial_r_alpha_x(colB[i][1], H.size(), p);
        int64_t evalB = evaluatePolynomial(buff_nB, beta1, p);
        vector<int64_t> buffB = calculatePolynomial_r_alpha_x(rowB[i][1], H.size(), p);
        evalB = (evalB *valB[i][1]) % p;
        if (evalB < 0) {
            evalB += p;
        }
        buffB = multiplyPolynomialByNumber(buffB, evalB, p);

        if(i > 0) {
            B_hat_M_hat = addPolynomials(B_hat_M_hat, buffB, p);
        }
        else {
            B_hat_M_hat = buffB;
        }
    }
    for (uint64_t i = 0; i < nonZeroRowsC.size(); i++) {
        vector<int64_t> buff_nC = calculatePolynomial_r_alpha_x(colC[i][1], H.size(), p);
        int64_t evalC = evaluatePolynomial(buff_nC, beta1, p);
        vector<int64_t> buffC = calculatePolynomial_r_alpha_x(rowC[i][1], H.size(), p);
        evalC = (evalC *valC[i][1]) % p;
        if (evalC < 0) {
            evalC += p;
        }
        buffC = multiplyPolynomialByNumber(buffC, evalC, p);

        if(i > 0) {
            C_hat_M_hat = addPolynomials(C_hat_M_hat, buffC, p);
        }
        else {
            C_hat_M_hat = buffC;
        }
    }
    storePolynomial(A_hat_M_hat, "A_hat_M_hat");
    storePolynomial(B_hat_M_hat, "B_hat_M_hat");
    storePolynomial(C_hat_M_hat, "C_hat_M_hat");

    
    vector<int64_t> eta_A_hat_M_hat = multiplyPolynomialByNumber(A_hat_M_hat, etaA, p);
    vector<int64_t> eta_B_hat_M_hat = multiplyPolynomialByNumber(B_hat_M_hat, etaB, p);
    vector<int64_t> eta_C_hat_M_hat = multiplyPolynomialByNumber(C_hat_M_hat, etaC, p);
    
    PrintPolynomial(eta_A_hat_M_hat, "eta_A_hat_M_hat: ");
    PrintPolynomial(eta_B_hat_M_hat, "eta_B_hat_M_hat: ");
    PrintPolynomial(eta_C_hat_M_hat, "eta_C_hat_M_hat: ");

    vector<int64_t> r_Sum_M_eta_M_M_hat_x_beta1 = multiplyPolynomials(addPolynomials(addPolynomials(eta_A_hat_M_hat, eta_B_hat_M_hat, p), eta_C_hat_M_hat, p), r_alpha_x, p);
    storePolynomial(r_Sum_M_eta_M_M_hat_x_beta1, "r_Sum_M_eta_M_M_hat_x_beta1");

    
    vector<int64_t> h_2_x = dividePolynomials(r_Sum_M_eta_M_M_hat_x_beta1, vH_x, p)[0];
    storePolynomial(h_2_x, "h2(x)");

    vector<int64_t> g_2_x = dividePolynomials(r_Sum_M_eta_M_M_hat_x_beta1, vH_x, p)[1];
    g_2_x.erase(g_2_x.begin());
    storePolynomial(g_2_x, "g2(x)");


    // int64_t beta2 = generateRandomNumber(H, p);
    int64_t beta2 = 80;
    cout << "beta2 = " << beta2 << endl;


    x_values.clear();
    y_values.clear();
    for (const auto& row : rowA) {
        if (!row.empty()) {
            x_values.push_back(row[0]);
            y_values.push_back(row[1]);
        }
    }
    setupLagrangePolynomial(x_values, y_values, p, "rowA(x)");
    
    x_values.clear();
    y_values.clear();
    for (const auto& col : colA) {
        if (!col.empty()) {
            x_values.push_back(col[0]);
            y_values.push_back(col[1]);
        }
    }
    setupLagrangePolynomial(x_values, y_values, p, "colA(x)");
    
    x_values.clear();
    y_values.clear();
    for (const auto& val : valA) {
        if (!val.empty()) {
            x_values.push_back(val[0]);
            y_values.push_back(val[1]);
        }
    }
    setupLagrangePolynomial(x_values, y_values, p, "valA(x)");
    
    x_values.clear();
    y_values.clear();
    for (const auto& row : rowB) {
        if (!row.empty()) {
            x_values.push_back(row[0]);
            y_values.push_back(row[1]);
        }
    }
    setupLagrangePolynomial(x_values, y_values, p, "rowB(x)");
    
    x_values.clear();
    y_values.clear();
    for (const auto& col : colB) {
        if (!col.empty()) {
            x_values.push_back(col[0]);
            y_values.push_back(col[1]);
        }
    }
    setupLagrangePolynomial(x_values, y_values, p, "colB(x)");
    
    x_values.clear();
    y_values.clear();
    for (const auto& val : valB) {
        if (!val.empty()) {
            x_values.push_back(val[0]);
            y_values.push_back(val[1]);
        }
    }
    setupLagrangePolynomial(x_values, y_values, p, "valB(x)");
    
    x_values.clear();
    y_values.clear();
    for (const auto& row : rowC) {
        if (!row.empty()) {
            x_values.push_back(row[0]);
            y_values.push_back(row[1]);
        }
    }
    setupLagrangePolynomial(x_values, y_values, p, "rowC(x)");
    
    x_values.clear();
    y_values.clear();
    for (const auto& col : colC) {
        if (!col.empty()) {
            x_values.push_back(col[0]);
            y_values.push_back(col[1]);
        }
    }
    setupLagrangePolynomial(x_values, y_values, p, "colC(x)");
    
    x_values.clear();
    y_values.clear();
    for (const auto& val : valC) {
        if (!val.empty()) {
            x_values.push_back(val[0]);
            y_values.push_back(val[1]);
        }
    }
    setupLagrangePolynomial(x_values, y_values, p, "valC(x)");

    int64_t vH_beta1 = evaluatePolynomial(vH_x, beta1, p);
    cout << "vH(beta1) = " << vH_beta1 << endl;

    int64_t vH_beta2 = evaluatePolynomial(vH_x, beta2, p);
    cout << "vH(beta2) = " << vH_beta2 << endl;





    vector<int64_t> rowA_x = parsePolynomial(LagrangeOutput["rowA(x)"]);
    vector<int64_t> colA_x = parsePolynomial(LagrangeOutput["colA(x)"]);
    vector<int64_t> valA_x = parsePolynomial(LagrangeOutput["valA(x)"]);
    
    vector<int64_t> rowB_x = parsePolynomial(LagrangeOutput["rowB(x)"]);
    vector<int64_t> colB_x = parsePolynomial(LagrangeOutput["colB(x)"]);
    vector<int64_t> valB_x = parsePolynomial(LagrangeOutput["valB(x)"]);
    
    vector<int64_t> rowC_x = parsePolynomial(LagrangeOutput["rowC(x)"]);
    vector<int64_t> colC_x = parsePolynomial(LagrangeOutput["colC(x)"]);
    vector<int64_t> valC_x = parsePolynomial(LagrangeOutput["valC(x)"]);

    ordered_json commitment ;
    commitment.clear();
    commitment["m"] = m;
    commitment["t"] = t;
    commitment["n"] = n;
    commitment["ComRowA"] = rowA_x;
    commitment["ComColA"] = colA_x;
    commitment["ComValA"] = valA_x;

    commitment["ComRowB"] = rowB_x;
    commitment["ComColB"] = colB_x;
    commitment["ComValB"] = valB_x;

    commitment["ComRowC"] = rowC_x;
    commitment["ComColC"] = colC_x;
    commitment["ComValC"] = valC_x;

    // Serialize JSON object to a string
    std::string commitmentString = commitment.dump();
    // Write JSON object to a file
    std::ofstream commitmentFile("commitment.json");
    if (commitmentFile.is_open()) {
        commitmentFile << commitmentString;
        commitmentFile.close();
        std::cout << "JSON data has been written to commitment.json\n";
    } else {
        std::cerr << "Error opening file for writing\n";
    }

    vector<int64_t> points_f_3 (K.size(), 0);
    int64_t sigma3 = 0;
    for (uint64_t i = 0; i < K.size(); i++) {
        int64_t sig3_A = 0;
        int64_t sig3_B = 0;
        int64_t sig3_C = 0;

        int64_t deA = (beta2 - evaluatePolynomial(rowA_x, K[i], p)) * (beta1 - evaluatePolynomial(colA_x, K[i], p)) % p;
        int64_t deB = (beta2 - evaluatePolynomial(rowB_x, K[i], p)) * (beta1 - evaluatePolynomial(colB_x, K[i], p)) % p;
        int64_t deC = (beta2 - evaluatePolynomial(rowC_x, K[i], p)) * (beta1 - evaluatePolynomial(colC_x, K[i], p)) % p;
        if (deA < 0) deA += p;
        if (deB < 0) deB += p;
        if (deC < 0) deC += p;

        sig3_A = etaA * ((vH_beta2 * vH_beta1 * evaluatePolynomial(valA_x, K[i], p)) * (mod_inverse(deA, p))) % p;
        sig3_B = etaB * ((vH_beta2 * vH_beta1 * evaluatePolynomial(valB_x, K[i], p)) * (mod_inverse(deB, p))) % p;
        sig3_C = etaC * ((vH_beta2 * vH_beta1 * evaluatePolynomial(valC_x, K[i], p)) * (mod_inverse(deC, p))) % p;
        
        points_f_3[i] = (sig3_A + sig3_B + sig3_C) % p;
        sigma3 += points_f_3[i];
        sigma3 %= p;
    }
    
    cout << "sigma3 = " << sigma3 << endl;
    
    vector<int64_t> poly_beta1 = {beta1};
    vector<int64_t> poly_beta2 = {beta2};
    
    vector<int64_t> poly_pi_a = multiplyPolynomials(subtractPolynomials(rowA_x, poly_beta2, p), subtractPolynomials(colA_x, poly_beta1, p), p);
    vector<int64_t> poly_pi_b = multiplyPolynomials(subtractPolynomials(rowB_x, poly_beta2, p), subtractPolynomials(colB_x, poly_beta1, p), p);
    vector<int64_t> poly_pi_c = multiplyPolynomials(subtractPolynomials(rowC_x, poly_beta2, p), subtractPolynomials(colC_x, poly_beta1, p), p);

    vector<int64_t> poly_etaA_vH_B2_vH_B1 = {(etaA * vH_beta2 * vH_beta1) % p};
    vector<int64_t> poly_etaB_vH_B2_vH_B1 = {(etaB * vH_beta2 * vH_beta1) % p};
    vector<int64_t> poly_etaC_vH_B2_vH_B1 = {(etaC * vH_beta2 * vH_beta1) % p};

    vector<int64_t> poly_sig_a = multiplyPolynomials(poly_etaA_vH_B2_vH_B1, valA_x, p);
    vector<int64_t> poly_sig_b = multiplyPolynomials(poly_etaB_vH_B2_vH_B1, valB_x, p);
    vector<int64_t> poly_sig_c = multiplyPolynomials(poly_etaC_vH_B2_vH_B1, valC_x, p);

    vector<int64_t> a_x = addPolynomials(addPolynomials(multiplyPolynomials(poly_sig_a, multiplyPolynomials(poly_pi_b, poly_pi_c, p), p), multiplyPolynomials(poly_sig_b, multiplyPolynomials(poly_pi_a, poly_pi_c, p), p), p), multiplyPolynomials(poly_sig_c, multiplyPolynomials(poly_pi_a, poly_pi_b, p), p), p);
    storePolynomial(a_x, "a(x)");
    
    vector<int64_t> b_x = multiplyPolynomials(multiplyPolynomials(poly_pi_a, poly_pi_b, p), poly_pi_c, p);
    storePolynomial(b_x, "b(x)");



    vector<int64_t> vK_x = createLinearPolynomial(K[0]);
    // Multiply (x - K) for all other K
    for (size_t i = 1; i < K.size(); i++) {
        vector<int64_t> nextPoly = createLinearPolynomial(K[i]);
        vK_x = multiplyPolynomials(vK_x, nextPoly, p);
    }
    storePolynomial(vK_x, "vK(x)");

    setupLagrangePolynomial(K, points_f_3, p, "poly_f_3(x)");
    vector<int64_t> poly_f_3x = parsePolynomial(LagrangeOutput["poly_f_3(x)"]);

    vector<int64_t> g_3_x = poly_f_3x;
    g_3_x.erase(g_3_x.begin());
    storePolynomial(g_3_x, "g3(x)");
    
    vector<int64_t> sigma_3_set_k;
    sigma_3_set_k.push_back((sigma3 * mod_inverse(K.size(), p)) % p);
    cout << "sigma_3_set_k: " << sigma_3_set_k[0] << endl;
    
    vector<int64_t> poly_f_3x_new = subtractPolynomials(poly_f_3x, sigma_3_set_k, p);
    storePolynomial(poly_f_3x_new, "f3(x)new");

    vector<int64_t> h_3_x = dividePolynomials(subtractPolynomials(a_x, multiplyPolynomials(b_x, addPolynomials(poly_f_3x_new, sigma_3_set_k, p), p), p), vK_x, p)[0];
    storePolynomial(h_3_x, "h3(x)");

    ordered_json proof;
    proof.clear();
    proof["P1AHP"] = sigma1;
    proof["P2AHP"] = w_hat_x;
    proof["P3AHP"] = z_hatA;
    proof["P4AHP"] = z_hatB;
    proof["P5AHP"] = z_hatC;
    proof["P6AHP"] = h_0_x;
    proof["P7AHP"] = s_x;
    proof["P8AHP"] = g_1_x;
    proof["P9AHP"] = h_1_x;
    proof["P10AHP"] = sigma2;
    proof["P11AHP"] = g_2_x;
    proof["P12AHP"] = h_2_x;
    proof["P13AHP"] = sigma3;
    proof["P14AHP"] = g_3_x;
    proof["P15AHP"] = h_3_x;
    proof["curve"] = "bn128";
    proof["protocol"] = "fidesv1";

    // Serialize JSON object to a string
    std::string proofString = proof.dump();

    // Use regex to format arrays correctly
    // proofString = std::regex_replace(proofString, std::regex("\\[\n\\s*([^\\]]*)\\s*\n\\]"), "[$1]");

    // Write JSON object to a file
    std::ofstream proofFile("proof.json");
    if (proofFile.is_open()) {
        proofFile << proofString;
        proofFile.close();
        std::cout << "JSON data has been written to proof.json\n";
    } else {
        std::cerr << "Error opening file for writing\n";
    }


    int64_t beta_3 = 5;



    // int64_t value = 1234567890123456789; // Example int64_t value

    // // Convert int64_t to byte array
    // unsigned char byteArray[sizeof(int64_t)];
    // int64ToByteArray(value, byteArray);

    // // Compute SHA256 hash
    // unsigned char hash[SHA256_DIGEST_LENGTH];
    // SHA256(byteArray, sizeof(byteArray), hash);

    // // Print the hash in hexadecimal format
    // std::cout << "SHA256 hash: ";
    // for (int i = 0; i < SHA256_DIGEST_LENGTH; ++i) {
    //     std::cout << std::hex << std::setw(2) << std::setfill('0') << static_cast<int>(hash[i]);
    // }
    // std::cout << std::endl;


}

void ZKP::AHPverify(int64_t mod) {
    int64_t alpha = 10;
    int64_t beta1 = 22;
    int64_t beta2 = 80;
    int64_t beta3 = 5;


    // Create a JSON object
    nlohmann::json jsonData;
    // Open the JSON file for reading
    std::ifstream proofFile("proof.json");
    if (!proofFile.is_open()) {
        std::cerr << "Could not open the file!" << std::endl;
    }

    // Parse the JSON data from the file
    proofFile >> jsonData;
    proofFile.close();

    // Store the data in appropriate types
    int64_t sigma1 = jsonData["P1AHP"].get<int64_t>();
    vector<int64_t> w_hat_x = jsonData["P2AHP"].get<vector<int64_t>>();
    vector<int64_t> z_hatA = jsonData["P3AHP"].get<vector<int64_t>>();
    vector<int64_t> z_hatB = jsonData["P4AHP"].get<vector<int64_t>>();
    vector<int64_t> z_hatC = jsonData["P5AHP"].get<vector<int64_t>>();
    vector<int64_t> h_0_x = jsonData["P6AHP"].get<vector<int64_t>>();
    vector<int64_t> s_x = jsonData["P7AHP"].get<vector<int64_t>>();
    vector<int64_t> g_1_x = jsonData["P8AHP"].get<vector<int64_t>>();
    vector<int64_t> h_1_x = jsonData["P9AHP"].get<vector<int64_t>>();
    int64_t sigma2 = jsonData["P10AHP"].get<int64_t>();
    vector<int64_t> g_2_x = jsonData["P11AHP"].get<vector<int64_t>>();
    vector<int64_t> h_2_x = jsonData["P12AHP"].get<vector<int64_t>>();
    int64_t sigma3 = jsonData["P13AHP"].get<int64_t>();
    vector<int64_t> g_3_x = jsonData["P14AHP"].get<vector<int64_t>>();
    vector<int64_t> h_3_x = jsonData["P15AHP"].get<vector<int64_t>>();
    string curve = jsonData["curve"];
    string protocol = jsonData["protocol"];


    std::ifstream commitmentFile("commitment.json");
    if (!commitmentFile.is_open()) {
        std::cerr << "Could not open the file!" << std::endl;
    }
    commitmentFile >> jsonData;
    commitmentFile.close();
    int64_t m = jsonData["m"].get<int64_t>();
    int64_t t = jsonData["t"].get<int64_t>();
    int64_t n = jsonData["n"].get<int64_t>();
    vector<int64_t> rowA_x = jsonData["ComRowA"].get<vector<int64_t>>();
    vector<int64_t> colA_x = jsonData["ComColA"].get<vector<int64_t>>();
    vector<int64_t> valA_x = jsonData["ComValA"].get<vector<int64_t>>();
    vector<int64_t> rowB_x = jsonData["ComRowB"].get<vector<int64_t>>();
    vector<int64_t> colB_x = jsonData["ComColB"].get<vector<int64_t>>();
    vector<int64_t> valB_x = jsonData["ComValB"].get<vector<int64_t>>();
    vector<int64_t> rowC_x = jsonData["ComRowC"].get<vector<int64_t>>();
    vector<int64_t> colC_x = jsonData["ComColC"].get<vector<int64_t>>();
    vector<int64_t> valC_x = jsonData["ComValC"].get<vector<int64_t>>();

    vector<int64_t> vK_x = parsePolynomial(LagrangeOutput["vK(x)"]);
    vector<int64_t> vH_x = parsePolynomial(LagrangeOutput["vH(x)"]);
    vector<int64_t> a_x = parsePolynomial(LagrangeOutput["a(x)"]);
    vector<int64_t> b_x = parsePolynomial(LagrangeOutput["b(x)"]);

    // vector<int64_t> r_alpha_x = parsePolynomial(LagrangeOutput["r(alpha, x)"]);
    vector<int64_t> r_alpha_x = calculatePolynomial_r_alpha_x(alpha, n, mod);

    vector<int64_t> Sum_M_eta_M_z_hat_M_x = parsePolynomial(LagrangeOutput["Sum_M_z_hatM(x)"]);
    vector<int64_t> z_hat_x = parsePolynomial(LagrangeOutput["z_hat(x)"]);


    

    bool verify = false;
    int64_t eq11 = (evaluatePolynomial(h_3_x, beta3, mod) * evaluatePolynomial(vK_x, beta3, mod)) % mod; // 22780
    int64_t eq12 = (evaluatePolynomial(a_x, beta3, mod) - ((evaluatePolynomial(b_x, beta3, mod) * (beta3 * evaluatePolynomial(g_3_x, beta3, mod) + (sigma3 * mod_inverse(m, mod))))) % mod); // 17 - 32 * (5*84 + 84/9)

    int64_t eq21 = (evaluatePolynomial(r_alpha_x, beta2, mod) * sigma3) % mod;
    int64_t eq22 = (evaluatePolynomial(h_2_x, beta2, mod) * evaluatePolynomial(vH_x, beta2, mod) + beta2 * evaluatePolynomial(g_2_x, beta2, mod) + sigma2 * mod_inverse(n, mod)) % mod;
    if (eq12 < 0) eq12 += mod;
    eq12 %= mod;

    int64_t eq31 = (evaluatePolynomial(s_x, beta1, mod) + evaluatePolynomial(r_alpha_x, beta1, mod) * evaluatePolynomial(Sum_M_eta_M_z_hat_M_x, beta1, mod) - sigma2 * evaluatePolynomial(z_hat_x, beta1, mod));
    int64_t eq32 = (evaluatePolynomial(h_1_x, beta1, mod) * evaluatePolynomial(vH_x, beta1, mod) + beta1 * evaluatePolynomial(g_1_x, beta1, mod) + sigma1 * mod_inverse(n, mod)) % mod;
    if (eq31 < 0) eq31 += mod;
    eq31 %= mod;

    int64_t eq41 = (evaluatePolynomial(z_hatA, beta1, mod) * evaluatePolynomial(z_hatB, beta1, mod) - evaluatePolynomial(z_hatC, beta1, mod));
    int64_t eq42 = (evaluatePolynomial(h_0_x, beta1, mod) * evaluatePolynomial(vH_x, beta1, mod)) % mod;
    if (eq41 < 0) eq41 += mod;
    eq41 %= mod;

    if (eq11 == eq12 && eq21 == eq22 && eq31 == eq32 && eq41 == eq42) {
        verify = true;
    }
    cout << "verify: " << verify << endl;
}

