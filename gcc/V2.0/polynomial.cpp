// Copyright 2024 FidesInnova.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.


#include "polynomial.h"
#include <iostream>
#include <unordered_map>
#include <random>
// #include <openssl/evp.h>
#include <iomanip>


#include <cstdint>
#include <cstring>
#include <vector>
#include <cinttypes>
#include <set>

#include <complex>
#include <cmath>
#include <algorithm>

int64_t Polynomial::power(int64_t base, int64_t exponent, int64_t p) {
  uint64_t result = 1;
  base = base % p;  // Handle base larger than p

  while (exponent > 0) {
    // If exponent is odd, multiply the result by base
    if (exponent % 2 == 1) {
      result = (result * base) % p;
    }
    // Square the base and reduce exponent by half
    base = (base * base) % p;
    exponent /= 2;
  }

  return result;
}

//These two functions Polynomial::power and Polynomial::pExp are the same, should delete one to be clear

// Function to compute the p exponentiation (a^b) % p
// int64_t Polynomial::pExp(int64_t a, int64_t b, int64_t p) {
//   int64_t result = 1;
//   while (b > 0) {
//     if (b & 1) result = (result * a) % p;
//     a = (a * a) % p;
//     b >>= 1;
//   }
//   if(result < 0) {
//     result += p;
//   }
//   return result;
// }
// Function to compute the p exponentiation (a^b) % p
int64_t Polynomial::pExp(int64_t a, int64_t b, int64_t p) {
  int64_t result = 1;
  a = a % p;
  while (b > 0) {
    if (b % 2 == 1) {  // If b is odd, multiply a with the result
      result = (result * a) % p;
    }
    b = b >> 1;         // Divide b by 2
    a = (a * a) % p;  // Change a to a^2
  }
  return result;
}

// Function to compute the p inverse using Fermat's Little Theorem
int64_t Polynomial::pInverse(int64_t a, int64_t p) {
  return pExp(a, p - 2, p);
}

// int64_t Polynomial::generateRandomNumber(const vector<int64_t>& H, int64_t p) {
//   int64_t randomNumber;
//   bool found;
//   do {
//     randomNumber = random(0, p);  // Generate a random number between 0 and p-1
//     found = false;
//     // Check if the generated number is in array H
//     for (size_t i = 0; i < H.size(); i++) {
//       if (H[i] == randomNumber) {
//         found = true;
//         break;
//       }
//     }
//   } while (found);
//   return randomNumber;
// }

int64_t Polynomial::generateRandomNumber(const std::vector<int64_t>& H, int64_t mod) {
    std::mt19937_64 rng(std::random_device{}());  // Use random_device to seed the generator
    std::uniform_int_distribution<int64_t> dist(0, mod - 1);
    
    int64_t randomNumber;
    do {
        randomNumber = dist(rng);
    } while (std::find(H.begin(), H.end(), randomNumber) != H.end());
    return randomNumber;
}

// Function to generate a random polynomial
vector<int64_t> Polynomial::generateRandomPolynomial(size_t numTerms, size_t maxDegree, int64_t p) {
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
        polynomial[index] = dis(gen) % p;
        if (polynomial[index] < 0) polynomial[index] += p;
    }

    return polynomial;
}

// Add two polynomials with p arithmetic
vector<int64_t> Polynomial::addPolynomials(const vector<int64_t>& poly1, const vector<int64_t>& poly2, int64_t p) {
  // Determine the size of the result polynomial (the max size of the two input polynomials)
  size_t maxSize = max(poly1.size(), poly2.size());
  vector<int64_t> result(maxSize, 0);

  // Add the polynomials
  for (size_t i = 0; i < maxSize; ++i) {
    if (i < poly1.size()) {
      result[i] = (result[i] + poly1[i]) % p;
    }
    if (i < poly2.size()) {
      result[i] = (result[i] + poly2[i]) % p;
    }
    // Ensure the result is non-negative
    if (result[i] < 0) {
      result[i] += p;
    }
  }
  return result;


    // size_t max_size = max(poly1.size(), poly2.size());
    // vector<int64_t> result(max_size, 0);

    // for (size_t i = 0; i < max_size; i++) {
    //     int64_t coeff1 = (i < poly1.size()) ? poly1[i] : 0;
    //     int64_t coeff2 = (i < poly2.size()) ? poly2[i] : 0;
    //     result[i] = (coeff1 + coeff2) % p;
    //     if (result[i] < 0) result[i] += p;  // Ensure positive modulo
    // }
    // return result;
}

// Subtract two polynomials with p arithmetic
vector<int64_t> Polynomial::subtractPolynomials(const vector<int64_t>& poly1, const vector<int64_t>& poly2, int64_t p) {
  size_t maxSize = max(poly1.size(), poly2.size());
  vector<int64_t> buf1 = poly1;
  vector<int64_t> buf2 = poly2;
  vector<int64_t> result(maxSize, 0);

  // Subtract the polynomials
  for (size_t i = 0; i < maxSize; i++) {
    int64_t val1 = (i < buf1.size()) ? buf1[i] : 0;
    int64_t val2 = (i < buf2.size()) ? buf2[i] : 0;

    result[i] = (val1 - val2) % p;
    if (result[i] < 0) {
      result[i] += p;  // Handle negative results
    }
  }

  return result;
}

// Function to multiply two polynomials

// Perform NTT or inverse NTT
void NTT(vector<int64_t>& a, bool invert, int64_t p, int64_t root) {
    int n = a.size();
    int64_t root_inv = Polynomial::pExp(root, p - 2, p);  // Inverse of root mod p
    int64_t root_pw = Polynomial::pExp(root, (p - 1) / n, p);  // root^(n/pw) mod p
    if (invert) root_pw = Polynomial::pExp(root_pw, p - 2, p);

    // Bit-reversal permutation
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        while (j & bit) {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if (i < j) swap(a[i], a[j]);
    }

    // NTT computation
    for (int len = 2; len <= n; len <<= 1) {
        int64_t wlen = Polynomial::pExp(root_pw, n / len, p);
        for (int i = 0; i < n; i += len) {
            int64_t w = 1;
            for (int j = 0; j < len / 2; j++) {
                int64_t u = a[i + j];
                int64_t v = (a[i + j + len / 2] * w) % p;
                a[i + j] = (u + v) % p;
                a[i + j + len / 2] = (u - v + p) % p;
                w = (w * wlen) % p;
            }
        }
    }

    // Scale for inverse NTT
    if (invert) {
        int64_t inv_n = Polynomial::pExp(n, p - 2, p);
        for (int64_t& x : a) x = (x * inv_n) % p;
    }
}

// Multiply two polynomials using NTT
// vector<int64_t> Polynomial::multiplyPolynomials(const vector<int64_t>& poly1, const vector<int64_t>& poly2, int64_t p) {
//     int64_t root = 15;
//     int n = 1;
//     while (n < poly1.size() + poly2.size() - 1) n <<= 1;

//     vector<int64_t> a(poly1.begin(), poly1.end());
//     vector<int64_t> b(poly2.begin(), poly2.end());
//     a.resize(n, 0);
//     b.resize(n, 0);

//     NTT(a, false, p, root);
//     NTT(b, false, p, root);

//     vector<int64_t> c(n);
//     for (int i = 0; i < n; i++) {
//         c[i] = (a[i] * b[i]) % p;
//     }

//     NTT(c, true, p, root);

//     // Remove trailing zeros
//     while (!c.empty() && c.back() == 0) c.pop_back();
//     return c;
// }


vector<int64_t> Polynomial::multiplyPolynomials(const vector<int64_t>& poly1, const vector<int64_t>& poly2, int64_t p) {
  vector<int64_t> result(poly1.size() + poly2.size() - 1, 0);

  for (size_t i = 0; i < poly1.size(); i++) {
    for (size_t j = 0; j < poly2.size(); j++) {
      result[i + j] = (result[i + j] + poly1[i] * poly2[j]) % p;
      if (result[i + j] < 0) {
        result[i + j] += p;
      }
    }
  }
  return result;
}


// Function to divide two polynomials
vector<vector<int64_t>> Polynomial::dividePolynomials(const vector<int64_t>& dividend, const vector<int64_t>& divisor, int64_t p) {
  vector<int64_t> quotient(dividend.size(), 0);  // Initialize with size equal to dividend size
  vector<int64_t> remainder = dividend;
  vector<vector<int64_t>> result;

  int64_t n = remainder.size();
  int64_t m = divisor.size();

  // If the divisor is larger than the dividend, the quotient is zero
  if (m > n) {
    quotient.resize(1, 0);  // Just one term, 0
    result.push_back(quotient);
    result.push_back(remainder);
    return result;
  }

  // Main loop for polynomial division
  while (remainder.size() >= divisor.size()) {
    // Degree of the current term in the quotient
    int64_t degreeDiff = remainder.size() - divisor.size();
    int64_t coef = (remainder.back() * pInverse(divisor.back(), p)) % p;

    // Add coefficient to the correct degree in the quotient
    quotient[degreeDiff] = coef;

    vector<int64_t> term(degreeDiff + 1, 0);
    term[degreeDiff] = coef;

    // Multiply the divisor by the current term in the quotient
    vector<int64_t> termTimesDivisor = multiplyPolynomials(term, divisor, p);
    termTimesDivisor.resize(remainder.size(), 0);  // Resize to match remainder size

    // Subtract to update the remainder
    remainder = subtractPolynomials(remainder, termTimesDivisor, p);

    // Remove leading zeros in the remainder
    while (remainder.size() > 0 && remainder.back() == 0) {
      remainder.pop_back();
    }
  }

  // Remove leading zeros in the quotient if necessary
  while (quotient.size() > 1 && quotient.back() == 0) {
    quotient.pop_back();
  }

  result.push_back(quotient);
  result.push_back(remainder);
  return result;
}

// Function to multiply a polynomial by a number
vector<int64_t> Polynomial::multiplyPolynomialByNumber(const vector<int64_t>& H, int64_t h, int64_t p) {
  vector<int64_t> result(H.size(), 0);  // Use long long to avoid overflow during multiplication

  for (int64_t i = 0; i < H.size(); i++) {
    result[i] = (H[i] * h) % p;

    // // If result becomes negative, convert it to a positive equivalent under pulo
    if (result[i] < 0) {
      result[i] += p;
    }
  }
  return result;
}

// Function to compute Lagrange basis polynomial L_i(x)
vector<int64_t> Polynomial::LagrangePolynomial(int64_t i, const vector<int64_t>& x_values, int64_t p) {
  int64_t n = x_values.size();
  vector<int64_t> result = {1};  // Start with 1 for the polynomial (constant term)
  vector<int64_t> inv_denominators(n);

  // Precompute modular inverses of (x_i - x_j) for all j ≠ i
  for (int64_t j = 0; j < n; j++) {
    if (j != i) {
      int64_t denominator = (x_values[i] + p - x_values[j]) % p;
      inv_denominators[j] = pInverse(denominator, p);
    }
  }

  for (int64_t j = 0; j < n; j++) {
    if (j != i) {
      vector<int64_t> term = { static_cast<int64_t>((p - x_values[j]) % p), 1 };  // (x - x_j)
      int64_t denominator = (x_values[i] + p - x_values[j]) % p;
      int64_t inv_denominator = inv_denominators[j];

      // Multiply the result by (x - x_j) / (x_i - x_j)
      vector<int64_t> temp = multiplyPolynomials(result, term, p);
      for (int64_t& coef : temp) {
        coef = (coef * inv_denominator) % p;
      }
      result = temp;
    }
  }
  return result;
}

// // Function to compute Lagrange polynomial(x, y)
// vector<int64_t> Polynomial::setupLagrangePolynomial(const vector<int64_t>& x_values, 
//                                                     const vector<int64_t>& y_values, 
//                                                     int64_t p, 
//                                                     const std::string& name) {
//     int64_t num_points = x_values.size();
//     vector<int64_t> polynomial(1, 0);  // Start with zero polynomial

//     for (int64_t i = 0; i < num_points; i++) {
//         if (y_values[i] != 0) {  // Process only non-zero y-values
//             vector<int64_t> Li = Polynomial::LagrangePolynomial(i, x_values, p);

//             // Scale Li by y_i modulo p
//             for (int64_t& coeff : Li) {
//                 coeff = (coeff * y_values[i]) % p;
//             }

//             // Add scaled Li to the final polynomial
//             polynomial = Polynomial::addPolynomials(polynomial, Li, p);
//         }
//     }

//     // Print the final polynomial
//     Polynomial::printPolynomial(polynomial, name);
//     return polynomial;
// }


vector<int64_t> Polynomial::setupLagrangePolynomial(const vector<int64_t>& x_values, const vector<int64_t>& y_values, int64_t p, const std::string& name) {
  // Automatically detect number of points
  int64_t num_points = x_values.size();

  // Compute polynomial(x) polynomial
  vector<int64_t> polynomial(1, 0);  // Start with a zero polynomial

  for (int64_t i = 0; i < num_points; i++) {
    if (y_values[i] != 0) {  // Only process non-zero y-values
      vector<int64_t> Li = LagrangePolynomial(i, x_values, p);
      // printPolynomial(Li, "L" + std::to_string(i + 1));

      // Multiply the L_i(x) by y_i and add to the final polynomial
      for (int64_t j = 0; j < Li.size(); j++) {
        if (j >= polynomial.size()) {
          polynomial.push_back(0);  // Ensure polynomial is large enough to accompate all terms
        }
        polynomial[j] = (polynomial[j] + y_values[i] * Li[j]) % p;
      }
    }
  }
  // Print the final polynomial(x) polynomial
  printPolynomial(polynomial, name);
  return polynomial;
}

// Function to parse the polynomial string and evaluate it
int64_t Polynomial::evaluatePolynomial(const vector<int64_t>& polynomial, int64_t x, int64_t p) {
  int64_t result = 0;

  for (size_t i = 0; i < polynomial.size(); i++) {

    int64_t termValue = polynomial[i] * power(x, i, p) % p;

    //handle when the number is negative
    if (termValue < 0) {
      termValue += p;
    }
    result += termValue;
    result %= p;
    if (result < 0) {
      result += p;
    }
  }
  return result;
}

// Function to compute the sum of polynomial evaluations at multiple points
int64_t Polynomial::sumOfEvaluations(const vector<int64_t>& poly, const vector<int64_t>& points, int64_t p) {
  //I add these lines of code to check if there is no elements in points
  // if (points.empty()) {
  //   cout << "Error: No points to evaluate" << endl;
  //   return 0;
  // }
  int64_t totalSum = 0;

  for (int64_t point : points) {
    int64_t evlpol = evaluatePolynomial(poly, point, p);
    totalSum = (totalSum + evlpol) % p;
  }
  totalSum %= p;

  if (totalSum < 0) totalSum += p;

  return totalSum;
}

// Function to create a polynomial for (x - root)
vector<int64_t> Polynomial::createLinearPolynomial(int64_t root) {
  return { -root, 1 };  // Represents (x - root)
}

// Function to calculate Polynomial r(α,x) = (alpha^n - x^n) / (alpha - x)
vector<int64_t> Polynomial::calculatePolynomial_r_alpha_x(int64_t alpha, int64_t n, int64_t p) {
  vector<int64_t> P(n, 0);

  // Calculate each term of the polynomial P(x)
  int64_t currentPowerOfAlpha = 1;  // alpha^0
  for (int64_t i = 0; i < n; ++i) {
    P[n - 1 - i] = currentPowerOfAlpha;  // alpha^(n-1-i)
    currentPowerOfAlpha = (currentPowerOfAlpha * alpha) % p;
  }

  return P;
}

// This function calculates the value of the polynomial at a specific point 𝑘 rather than returning the coefficients.
int64_t Polynomial::calculatePolynomial_r_alpha_k(int64_t alpha, int64_t k, int64_t n, int64_t p) {
  int64_t result = 1;
  result = (power(alpha, n, p) - power(k, n, p)) % p;
  int64_t buff = (alpha - k) % p;
  if (buff < 0) {
    buff += p;
  }
  result *= pInverse(buff, p);
  ;
  result %= p;
  return result;
}

// Function to expand polynomials given the roots
vector<int64_t> Polynomial::expandPolynomials(const vector<int64_t>& roots, int64_t p) {
  vector<int64_t> result = { 1 };  // Start with the polynomial "1"

  for (int64_t root : roots) {
    // Multiply the current result polynomial by (x - root)
    vector<int64_t> temp(result.size() + 1, 0);
    for (size_t i = 0; i < result.size(); i++) {
      temp[i] += result[i];  // x^n term
      temp[i] %= p;
      temp[i + 1] -= result[i] * root;  // -root * x^(n-1) term
      temp[i + 1] %= p;
      if (temp[i + 1] < 0) temp[i + 1] += p;
    }
    result = temp;
  }
  reverse(result.begin(), result.end());
  return result;
}

// Function to print polynomial in serial
void Polynomial::printPolynomial(const vector<int64_t>& coefficients, const std::string& name) {
  // Iterate through the coefficients and print each term of the polynomial
  cout << name  << " = ";
  bool first = true;
  for (int64_t i = coefficients.size() - 1; i >= 0; i--) {
    if (coefficients[i] == 0) continue;  // Skip zero coefficients

    // Print the sign for all terms except the first
    if (!first) {
      if (coefficients[i] > 0) {
        cout << " + ";
      } else {
        cout << " - ";
      }
    } else {
      first = false;
    }

    // Print the absolute value of the coefficient
    cout << abs(coefficients[i]);

    // Print the variable and the exponent
    cout << "x^" << i;
  }
  cout << endl;
}

// Utility functions for trimming and removing commas
std::string Polynomial::trim(const std::string& str) {
    size_t first = str.find_first_not_of(' ');
    if (first == std::string::npos)
        return "";
    size_t last = str.find_last_not_of(' ');
    return str.substr(first, last - first + 1);
}

// Helper function to Remove commas from a string
std::string Polynomial::removeCommas(const std::string& str) {
    size_t first = str.find_first_not_of(',');
    if (first == std::string::npos)
        return "";
    size_t last = str.find_last_not_of(',');
    return str.substr(first, last - first + 1);
}

void Polynomial::printMatrix(vector<vector<int64_t>>& matrix, const std::string& name) {
  cout << "Matrix " << name << ":" << endl;
  for (const auto& row : matrix) {
    for (int64_t val : row) {
      cout << val << " ";
    }
    cout << endl;
  }
}


// Function to get the row indices of non-zero entries in matrix
vector<vector<int64_t>> Polynomial::getNonZeroRows(const vector<vector<int64_t>>& matrix) {
  vector<vector<int64_t>> nonZeroRows(2);
  for (int64_t i = 0; i < matrix.size(); i++) {
    for (int64_t j = 0; j < matrix[i].size(); j++) {
      if (matrix[i][j] != 0) {
        nonZeroRows[0].push_back(i);  // Storing row index
        nonZeroRows[1].push_back(matrix[i][j]);
      }
    }
  }
  return nonZeroRows;
}

// Function to get the col indices of non-zero entries in matrix
vector<vector<int64_t>> Polynomial::getNonZeroCols(const vector<vector<int64_t>>& matrix) {
  vector<vector<int64_t>> nonZeroCols(2);
  for (int64_t i = 0; i < matrix.size(); i++) {
    for (int64_t j = 0; j < matrix[i].size(); j++) {
      if (matrix[i][j] != 0) {
        nonZeroCols[0].push_back(j);  // Storing col index
        nonZeroCols[1].push_back(matrix[i][j]);
      }
    }
  }
  return nonZeroCols;
}


// Function to create the mapping
vector<vector<int64_t>> Polynomial::createMapping(const vector<int64_t>& K, const vector<int64_t>& H, const vector<vector<int64_t>>& nonZero) {
  vector<vector<int64_t>> row(2);
  for (int64_t i = 0; i < nonZero[0].size(); i++) {
    row[0].push_back(K[i]);
    row[1].push_back(H[nonZero[0][i]]);
  }
  for (int64_t i = nonZero[0].size(); i < K.size(); i++) {
    row[0].push_back(K[i]);
    row[1].push_back(H[i % H.size()]);
  }

  return row;
}

// Function to print the mapping
void Polynomial::printMapping(vector<vector<int64_t>>& row, const std::string& name) {
  for (int64_t i = 0; i < row[0].size(); i++) {
    cout << name << "(" << row[0][i] << ") = " << row[1][i] << endl;
  }
  cout << endl;
}

// Function to create the val mapping
vector<vector<int64_t>> Polynomial::valMapping(const vector<int64_t>& K, const vector<int64_t>& H, vector<vector<int64_t>>& nonZeroRows, vector<vector<int64_t>>& nonZeroCols, int64_t p) {
  vector<vector<int64_t>> val(2);

  for (int64_t i = 0; i < K.size(); i++) {
    if (i < nonZeroRows[0].size()) {
      val[0].push_back(K[i]);
      val[1].push_back((nonZeroRows[1][i] * Polynomial::pInverse(((H.size() * Polynomial::power(H[nonZeroRows[0][i]], H.size() - 1, p)) % p) * ((H.size() * Polynomial::power(H[nonZeroCols[0][i]], H.size() - 1, p)) % p), p)) % p);
    } else {
      val[0].push_back(K[i]);
      val[1].push_back(0);
    }
  }

  return val;
}

// Function to calculate log in p
int64_t Polynomial::log_p(int64_t a, int64_t b, int64_t p) {
  int64_t m = static_cast<int64_t>(ceil(sqrt(static_cast<double>(p - 1))));

  // Precompute a^i (p P) for i = 0..m and store in a map
  std::unordered_map<int64_t, int64_t> tbl;
  for (int64_t i = 0; i < m; ++i) {
    tbl[Polynomial::pExp(a, i, p)] = i;
  }

  int64_t c = Polynomial::pExp(a, m * (p - 2), p);  // c = a^(-m) p P

  // Check if we can find the solution in the baby-step giant-step manner
  for (int64_t j = 0; j < m; ++j) {
    int64_t y = (b * Polynomial::pExp(c, j, p)) % p;
    if (tbl.find(y) != tbl.end()) {
      int64_t num = tbl[y];
      return j * m + num;
    }
  }
  return 0;  // Return 0 if no solution is found
}

// Function to calculate e_func in p
int64_t Polynomial::e_func(int64_t a, int64_t b, int64_t g, int64_t p) {
  int64_t buf1 = (a * Polynomial::pInverse(g, p)) % p;
  int64_t buf2 = (b * Polynomial::pInverse(g, p)) % p;
  // return (3 * (buf1 * buf2)) % p;
  return (buf1 * buf2) % p;
}

  // Function to calculate KZG in p
int64_t Polynomial::KZG_Commitment(vector<int64_t> a, vector<int64_t> b, int64_t p) {
  int64_t res = 0;
  for (int64_t i = 0; i < b.size(); i++) {
    res += (a[i] * b[i]) % p;
    res %= p;
  }
  return res;
}



// Function to compute the SHA-256 hash of an int64_t and return the lower 4 bytes as int64_t, applying a modulo operation
int64_t Polynomial::hashAndExtractLower4Bytes(int64_t inputNumber, int64_t p) {
  char* inputData = (char*)malloc(21);
  snprintf(inputData, 21, "%" PRId64, inputNumber);

  // Compute the hash
  string hashStr = Polynomial::SHA256(inputData);
  string last4BytesHex = hashStr.substr(56, 8);
  
  int64_t result = 0;
  std::stringstream ss;
  ss << std::hex << last4BytesHex;
  ss >> result;

  // Apply modulo operation
  result = result % p;

  // Ensure positive result for modulo
  if (result < 0) {
    result += p;
  }

  return result;
}





#define uchar unsigned char
#define uint unsigned int

#define DBL_INT_ADD(a,b,c) if (a > 0xffffffff - (c)) ++b; a += c;
#define ROTLEFT(a,b) (((a) << (b)) | ((a) >> (32-(b))))
#define ROTRIGHT(a,b) (((a) >> (b)) | ((a) << (32-(b))))

#define CH(x,y,z) (((x) & (y)) ^ (~(x) & (z)))
#define MAJ(x,y,z) (((x) & (y)) ^ ((x) & (z)) ^ ((y) & (z)))
#define EP0(x) (ROTRIGHT(x,2) ^ ROTRIGHT(x,13) ^ ROTRIGHT(x,22))
#define EP1(x) (ROTRIGHT(x,6) ^ ROTRIGHT(x,11) ^ ROTRIGHT(x,25))
#define SIG0(x) (ROTRIGHT(x,7) ^ ROTRIGHT(x,18) ^ ((x) >> 3))
#define SIG1(x) (ROTRIGHT(x,17) ^ ROTRIGHT(x,19) ^ ((x) >> 10))

typedef struct {
	uchar data[64];
	uint datalen;
	uint bitlen[2];
	uint state[8];
} SHA256_CTX;

uint k[64] = {
	0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
	0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
	0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
	0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
	0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
	0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
	0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
	0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2
};

void SHA256Transform(SHA256_CTX *ctx, uchar data[])
{
	uint a, b, c, d, e, f, g, h, i, j, t1, t2, m[64];

	for (i = 0, j = 0; i < 16; ++i, j += 4)
		m[i] = (data[j] << 24) | (data[j + 1] << 16) | (data[j + 2] << 8) | (data[j + 3]);
	for (; i < 64; ++i)
		m[i] = SIG1(m[i - 2]) + m[i - 7] + SIG0(m[i - 15]) + m[i - 16];

	a = ctx->state[0];
	b = ctx->state[1];
	c = ctx->state[2];
	d = ctx->state[3];
	e = ctx->state[4];
	f = ctx->state[5];
	g = ctx->state[6];
	h = ctx->state[7];

	for (i = 0; i < 64; ++i) {
		t1 = h + EP1(e) + CH(e, f, g) + k[i] + m[i];
		t2 = EP0(a) + MAJ(a, b, c);
		h = g;
		g = f;
		f = e;
		e = d + t1;
		d = c;
		c = b;
		b = a;
		a = t1 + t2;
	}

	ctx->state[0] += a;
	ctx->state[1] += b;
	ctx->state[2] += c;
	ctx->state[3] += d;
	ctx->state[4] += e;
	ctx->state[5] += f;
	ctx->state[6] += g;
	ctx->state[7] += h;
}

void SHA256Init(SHA256_CTX *ctx)
{
	ctx->datalen = 0;
	ctx->bitlen[0] = 0;
	ctx->bitlen[1] = 0;
	ctx->state[0] = 0x6a09e667;
	ctx->state[1] = 0xbb67ae85;
	ctx->state[2] = 0x3c6ef372;
	ctx->state[3] = 0xa54ff53a;
	ctx->state[4] = 0x510e527f;
	ctx->state[5] = 0x9b05688c;
	ctx->state[6] = 0x1f83d9ab;
	ctx->state[7] = 0x5be0cd19;
}

void SHA256Update(SHA256_CTX *ctx, uchar data[], uint len)
{
	for (uint i = 0; i < len; ++i) {
		ctx->data[ctx->datalen] = data[i];
		ctx->datalen++;
		if (ctx->datalen == 64) {
			SHA256Transform(ctx, ctx->data);
			DBL_INT_ADD(ctx->bitlen[0], ctx->bitlen[1], 512);
			ctx->datalen = 0;
		}
	}
}

void SHA256Final(SHA256_CTX *ctx, uchar hash[])
{
	uint i = ctx->datalen;

	if (ctx->datalen < 56) {
		ctx->data[i++] = 0x80;

		while (i < 56)
			ctx->data[i++] = 0x00;
	}
	else {
		ctx->data[i++] = 0x80;

		while (i < 64)
			ctx->data[i++] = 0x00;

		SHA256Transform(ctx, ctx->data);
		memset(ctx->data, 0, 56);
	}

	DBL_INT_ADD(ctx->bitlen[0], ctx->bitlen[1], ctx->datalen * 8);
	ctx->data[63] = ctx->bitlen[0];
	ctx->data[62] = ctx->bitlen[0] >> 8;
	ctx->data[61] = ctx->bitlen[0] >> 16;
	ctx->data[60] = ctx->bitlen[0] >> 24;
	ctx->data[59] = ctx->bitlen[1];
	ctx->data[58] = ctx->bitlen[1] >> 8;
	ctx->data[57] = ctx->bitlen[1] >> 16;
	ctx->data[56] = ctx->bitlen[1] >> 24;
	SHA256Transform(ctx, ctx->data);

	for (i = 0; i < 4; ++i) {
		hash[i] = (ctx->state[0] >> (24 - i * 8)) & 0x000000ff;
		hash[i + 4] = (ctx->state[1] >> (24 - i * 8)) & 0x000000ff;
		hash[i + 8] = (ctx->state[2] >> (24 - i * 8)) & 0x000000ff;
		hash[i + 12] = (ctx->state[3] >> (24 - i * 8)) & 0x000000ff;
		hash[i + 16] = (ctx->state[4] >> (24 - i * 8)) & 0x000000ff;
		hash[i + 20] = (ctx->state[5] >> (24 - i * 8)) & 0x000000ff;
		hash[i + 24] = (ctx->state[6] >> (24 - i * 8)) & 0x000000ff;
		hash[i + 28] = (ctx->state[7] >> (24 - i * 8)) & 0x000000ff;
	}
}

string Polynomial::SHA256(char* data) {
	int strLen = strlen(data);
	SHA256_CTX ctx;
	unsigned char hash[32];
	string hashStr = "";

	SHA256Init(&ctx);
	SHA256Update(&ctx, (unsigned char*)data, strLen);
	SHA256Final(&ctx, hash);

	char s[3];
	for (int i = 0; i < 32; i++) {
		sprintf(s, "%02x", hash[i]);
		hashStr += s;
	}

	return hashStr;
}
