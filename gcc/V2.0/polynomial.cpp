#include "polynomial.h"
#include <iostream>
#include <unordered_map>
#include <random>


int64_t Polynomial::power(int64_t base, int64_t exponent, int64_t p) {
  int64_t result = 1;
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

int64_t generateRandomNumber(const std::vector<int64_t>& H, int64_t mod) {
    std::mt19937_64 rng(std::random_device{}());  // Use random_device to seed the generator
    std::uniform_int_distribution<int64_t> dist(0, mod - 1);
    
    int64_t randomNumber;
    do {
        randomNumber = dist(rng);
    } while (std::find(H.begin(), H.end(), randomNumber) != H.end());

    return randomNumber;
}

// Add two polynomials with p arithmetic
vector<int64_t> Polynomial::addPolynomials(const vector<int64_t>& poly1, const vector<int64_t>& poly2, int64_t p) {
  // // I add this to check if p is positive for clarity usability of the function
  //   if (p <= 0) {
  //       cout << "Error: pulus must be a positive integer." << endl;
  //       return {};  // Return an empty vector to indicate an error
  //   }

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

//didn't test yet
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
  vector<int64_t> result = { 1 };  // Start with 1 for the polynomial (constant term)

  for (int64_t j = 0; j < n; j++) {
    if (j != i) {
      vector<int64_t> term = { static_cast<int64_t>((p - x_values[j]) % p), 1 };  // (x - x_j)
      int64_t denominator = (x_values[i] + p - x_values[j]) % p;

      // Efficient calculation of the p inverse
      int64_t inv_denominator = pInverse(denominator, p);

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

// Function to compute Lagrange polynomial(x, y)
vector<int64_t> Polynomial::setupLagrangePolynomial(const vector<int64_t> x_values, const vector<int64_t> y_values, int64_t p, const std::string& name) {
  // Automatically detect number of points
  int64_t num_points = x_values.size();

  // Compute polynomial(x) polynomial
  vector<int64_t> polynomial(1, 0);  // Start with a zero polynomial

  for (int64_t i = 0; i < num_points; i++) {
    if (y_values[i] != 0) {  // Only process non-zero y-values
      vector<int64_t> Li = LagrangePolynomial(i, x_values, p);
      printPolynomial(Li, "L" + std::to_string(i + 1));

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
  //add this to check if the p is 0, if it is return 0 which save unnecessary computation
  // if (p == 1) return 0;  // Early exit if p is 1

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
  return (3 * (buf1 * buf2)) % p;
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
