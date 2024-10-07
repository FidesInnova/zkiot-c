#include "polynomial.h"


int64_t Polynomial::power(int64_t base, int64_t exponent, int64_t mod) {
  int64_t result = 1;
  base = base % mod;  // Handle base larger than mod

  while (exponent > 0) {
    // If exponent is odd, multiply the result by base
    if (exponent % 2 == 1) {
      result = (result * base) % mod;
    }
    // Square the base and reduce exponent by half
    base = (base * base) % mod;
    exponent /= 2;
  }

  return result;
}

// Function to compute the modular exponentiation (a^b) % mod
int64_t Polynomial::modExp(int64_t a, int64_t b, int64_t mod) {
  int64_t result = 1;
  a = a % mod;
  while (b > 0) {
    if (b % 2 == 1) {  // If b is odd, multiply a with the result
      result = (result * a) % mod;
    }
    b = b >> 1;         // Divide b by 2
    a = (a * a) % mod;  // Change a to a^2
  }
  return result;
}

// Function to compute the modular inverse using Fermat's Little Theorem
int64_t Polynomial::modInverse(int64_t a, int64_t mod) {
  return modExp(a, mod - 2, mod);
}

int64_t Polynomial::generateRandomNumber(const vector<int64_t>& H, int64_t mod) {
  int64_t randomNumber;
  bool found;
  do {
    randomNumber = random(0, mod);  // Generate a random number between 0 and mod-1
    found = false;
    // Check if the generated number is in array H
    for (size_t i = 0; i < H.size(); i++) {
      if (H[i] == randomNumber) {
        found = true;
        break;
      }
    }
  } while (found);
  return randomNumber;
}

// Add two polynomials with modular arithmetic
vector<int64_t> Polynomial::addPolynomials(const vector<int64_t>& poly1, const vector<int64_t>& poly2, int64_t mod) {
  // Determine the size of the result polynomial (the max size of the two input polynomials)
  size_t maxSize = max(poly1.size(), poly2.size());
  vector<int64_t> result(maxSize, 0);

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

// Subtract two polynomials with modular arithmetic
vector<int64_t> Polynomial::subtractPolynomials(const vector<int64_t>& poly1, const vector<int64_t>& poly2, int64_t mod) {
  size_t maxSize = max(poly1.size(), poly2.size());
  vector<int64_t> buf1 = poly1;
  vector<int64_t> buf2 = poly2;
  vector<int64_t> result(maxSize, 0);

  // Subtract the polynomials
  for (size_t i = 0; i < maxSize; i++) {
    int64_t val1 = (i < buf1.size()) ? buf1[i] : 0;
    int64_t val2 = (i < buf2.size()) ? buf2[i] : 0;

    result[i] = (val1 - val2) % mod;
    if (result[i] < 0) {
      result[i] += mod;  // Handle negative results
    }
  }

  return result;
}

// Function to multiply two polynomials
vector<int64_t> Polynomial::multiplyPolynomials(const vector<int64_t>& poly1, const vector<int64_t>& poly2, int64_t mod) {
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

// Function to divide two polynomials
vector<vector<int64_t>> Polynomial::dividePolynomials(const vector<int64_t>& dividend, const vector<int64_t>& divisor, int64_t mod) {
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
    int64_t coef = (remainder.back() * modInverse(divisor.back(), mod)) % mod;

    // Add coefficient to the correct degree in the quotient
    quotient[degreeDiff] = coef;

    vector<int64_t> term(degreeDiff + 1, 0);
    term[degreeDiff] = coef;

    // Multiply the divisor by the current term in the quotient
    vector<int64_t> termTimesDivisor = multiplyPolynomials(term, divisor, mod);
    termTimesDivisor.resize(remainder.size(), 0);  // Resize to match remainder size

    // Subtract to update the remainder
    remainder = subtractPolynomials(remainder, termTimesDivisor, mod);

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
vector<int64_t> Polynomial::multiplyPolynomialByNumber(const vector<int64_t>& H, int64_t h, int64_t mod) {
  vector<int64_t> result(H.size(), 0);  // Use long long to avoid overflow during multiplication

  for (int64_t i = 0; i < H.size(); i++) {
    result[i] = (H[i] * h) % mod;

    // // If result becomes negative, convert it to a positive equivalent under modulo
    if (result[i] < 0) {
      result[i] += mod;
    }
  }
  return result;
}

// Function to compute Lagrange basis polynomial L_i(x)
vector<int64_t> Polynomial::LagrangePolynomial(int64_t i, const vector<int64_t>& x_values, int64_t mod) {
  int64_t n = x_values.size();
  vector<int64_t> result = { 1 };  // Start with 1 for the polynomial (constant term)

  for (int64_t j = 0; j < n; j++) {
    if (j != i) {
      vector<int64_t> term = { static_cast<int64_t>((mod - x_values[j]) % mod), 1 };  // (x - x_j)
      int64_t denominator = (x_values[i] + mod - x_values[j]) % mod;

      // Efficient calculation of the modular inverse
      int64_t inv_denominator = modInverse(denominator, mod);

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

// Function to compute Lagrange polynomial(x, y)
vector<int64_t> Polynomial::setupLagrangePolynomial(const vector<int64_t> x_values, const vector<int64_t> y_values, int64_t mod, const String& name) {
  // Automatically detect number of points
  int64_t num_points = x_values.size();

  // Compute polynomial(x) polynomial
  vector<int64_t> polynomial(1, 0);  // Start with a zero polynomial

  for (int64_t i = 0; i < num_points; i++) {
    if (y_values[i] != 0) {  // Only process non-zero y-values
      vector<int64_t> Li = LagrangePolynomial(i, x_values, mod);
      printPolynomial(Li, "L" + String(i + 1));

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
  printPolynomial(polynomial, name);
  return polynomial;
}

// Function to parse the polynomial string and evaluate it
int64_t Polynomial::evaluatePolynomial(const vector<int64_t>& polynomial, int64_t x, int64_t mod) {
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
int64_t Polynomial::sumOfEvaluations(const vector<int64_t>& poly, const vector<int64_t>& points, int64_t mod) {
  int64_t totalSum = 0;

  for (int64_t point : points) {
    int64_t evlpol = evaluatePolynomial(poly, point, mod);
    totalSum = (totalSum + evlpol) % mod;
  }
  totalSum %= mod;

  if (totalSum < 0) totalSum += mod;

  return totalSum;
}

// Function to create a polynomial for (x - root)
vector<int64_t> Polynomial::createLinearPolynomial(int64_t root) {
  return { -root, 1 };  // Represents (x - root)
}

// Function to calculate Polynomial r(α,x) = (alpha^n - x^n) / (alpha - x)
vector<int64_t> Polynomial::calculatePolynomial_r_alpha_x(int64_t alpha, int64_t n, int64_t mod) {
  vector<int64_t> P(n, 0);

  // Calculate each term of the polynomial P(x)
  int64_t currentPowerOfAlpha = 1;  // alpha^0
  for (int64_t i = 0; i < n; ++i) {
    P[n - 1 - i] = currentPowerOfAlpha;  // alpha^(n-1-i)
    currentPowerOfAlpha = (currentPowerOfAlpha * alpha) % mod;
  }

  return P;
}

// Function to calculate Polynomial r(α,x) = (alpha^n - x^n) / (alpha - x)
int64_t Polynomial::calculatePolynomial_r_alpha_k(int64_t alpha, int64_t k, int64_t n, int64_t mod) {
  int64_t result = 1;
  result = (power(alpha, n, mod) - power(k, n, mod)) % mod;
  int64_t buff = (alpha - k) % mod;
  if (buff < 0) {
    buff += mod;
  }
  result *= modInverse(buff, mod);
  ;
  result %= mod;
  return result;
}

// Function to expand polynomials given the roots
vector<int64_t> Polynomial::expandPolynomials(const vector<int64_t>& roots, int64_t mod) {
  vector<int64_t> result = { 1 };  // Start with the polynomial "1"

  for (int64_t root : roots) {
    // Multiply the current result polynomial by (x - root)
    vector<int64_t> temp(result.size() + 1, 0);
    for (size_t i = 0; i < result.size(); i++) {
      temp[i] += result[i];  // x^n term
      temp[i] %= mod;
      temp[i + 1] -= result[i] * root;  // -root * x^(n-1) term
      temp[i + 1] %= mod;
      if (temp[i + 1] < 0) temp[i + 1] += mod;
    }
    result = temp;
  }
  reverse(result.begin(), result.end());
  return result;
}

// Function to print polynomial in serial
void Polynomial::printPolynomial(const vector<int64_t>& coefficients, const String& name) {
  // Iterate through the coefficients and print each term of the polynomial
  Serial.print(name);
  Serial.print(" = ");
  bool first = true;
  for (int64_t i = coefficients.size() - 1; i >= 0; i--) {
    if (coefficients[i] == 0) continue;  // Skip zero coefficients

    // Print the sign for all terms except the first
    if (!first) {
      if (coefficients[i] > 0) {
        Serial.print(" + ");
      } else {
        Serial.print(" - ");
      }
    } else {
      first = false;
    }

    // Print the absolute value of the coefficient
    Serial.print(abs(coefficients[i]));

    // Print the variable and the exponent
    Serial.print("x^");
    Serial.print(i);
  }
  Serial.println("");
}

// Utility functions for trimming and removing commas
String Polynomial::trim(const String& str) {
  int64_t start = 0;
  while (start < str.length() && isspace(str[start])) start++;
  int64_t end = str.length() - 1;
  while (end >= 0 && isspace(str[end])) end--;
  return str.substring(start, end + 1);
}

String Polynomial::removeCommas(const String& str) {
  String result = str;
  result.replace(",", "");
  return result;
}

void Polynomial::printMatrix(vector<vector<int64_t>>& matrix, const String& name) {
  Serial.print("Matrix ");
  Serial.print(name);
  Serial.println(":");
  for (const auto& row : matrix) {
    for (int64_t val : row) {
      Serial.print(val);
      Serial.print(" ");
    }
    Serial.println("");
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
    // row[1].push_back(H[i % H.size()]);
  }

  return row;
}

// Function to print the mapping
void Polynomial::printMapping(vector<vector<int64_t>>& row, const String& name) {
  for (int64_t i = 0; i < row[0].size(); i++) {
    Serial.print(name);
    Serial.print("(");
    Serial.print(row[0][i]);
    Serial.print(") = ");
    Serial.println(row[1][i]);
  }
  Serial.println("");
}

// Function to create the val mapping
vector<vector<int64_t>> Polynomial::valMapping(const vector<int64_t>& K, const vector<int64_t>& H, vector<vector<int64_t>>& nonZeroRows, vector<vector<int64_t>>& nonZeroCols, int64_t mod) {
  vector<vector<int64_t>> val(2);

  for (int64_t i = 0; i < K.size(); i++) {
    if (i < nonZeroRows[0].size()) {
      val[0].push_back(K[i]);
      val[1].push_back((nonZeroRows[1][i] * Polynomial::modInverse(((H.size() * Polynomial::power(H[nonZeroRows[0][i]], H.size() - 1, mod)) % mod) * ((H.size() * Polynomial::power(H[nonZeroCols[0][i]], H.size() - 1, mod)) % mod), mod)) % mod);
    } else {
      val[0].push_back(K[i]);
      val[1].push_back(0);
    }
  }

  return val;
}

// Function to calculate log in mod
int64_t Polynomial::log_mod(int64_t a, int64_t b, int64_t mod) {
  int64_t m = static_cast<int64_t>(ceil(sqrt(static_cast<double>(mod - 1))));

  // Precompute a^i (mod P) for i = 0..m and store in a map
  std::unordered_map<int64_t, int64_t> tbl;
  for (int64_t i = 0; i < m; ++i) {
    tbl[Polynomial::modExp(a, i, mod)] = i;
  }

  int64_t c = Polynomial::modExp(a, m * (mod - 2), mod);  // c = a^(-m) mod P

  // Check if we can find the solution in the baby-step giant-step manner
  for (int64_t j = 0; j < m; ++j) {
    int64_t y = (b * Polynomial::modExp(c, j, mod)) % mod;
    if (tbl.find(y) != tbl.end()) {
      int64_t num = tbl[y];
      return j * m + num;
    }
  }

  return 0;  // Return 0 if no solution is found
}

// Function to calculate e_func in mod
int64_t Polynomial::e_func(int64_t a, int64_t b, int64_t g, int64_t mod) {
  int64_t buf1 = Polynomial::log_mod(g, a, mod);
  int64_t buf2 = Polynomial::log_mod(g, b, mod);
  Serial.print("buf1 = ");
  Serial.println(buf1);
  Serial.print("buf2 = ");
  Serial.println(buf2);
  return (Polynomial::power(3, (buf1 * buf2), mod));
}

