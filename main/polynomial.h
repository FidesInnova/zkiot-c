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

#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <Arduino.h>
#include <vector>
#include <cstdint>
#include <algorithm>

using namespace std;

class Polynomial {
public:
  // Function to compute the modular exponentiation (base^exponent) % mod
  static int64_t power(int64_t base, int64_t exponent, int64_t mod);

  // Function to compute the modular exponentiation (a^b) % mod
  static int64_t modExp(int64_t a, int64_t b, int64_t mod);

  // Function to compute the modular inverse using Fermat's Little Theorem
  static int64_t modInverse(int64_t a, int64_t mod);

  // Function to generate a random number in mod
  static int64_t generateRandomNumber(const vector<int64_t>& H, int64_t mod);

  // Add two polynomials with modular arithmetic
  static vector<int64_t> addPolynomials(const vector<int64_t>& poly1, const vector<int64_t>& poly2, int64_t mod);

  // Subtract two polynomials with modular arithmetic
  static vector<int64_t> subtractPolynomials(const vector<int64_t>& poly1, const vector<int64_t>& poly2, int64_t mod);

  // Function to multiply two polynomials
  static vector<int64_t> multiplyPolynomials(const vector<int64_t>& poly1, const vector<int64_t>& poly2, int64_t mod);

  // Function to divide two polynomials
  static vector<vector<int64_t>> dividePolynomials(const vector<int64_t>& dividend, const vector<int64_t>& divisor, int64_t mod);

  // Function to multiply a polynomial by a number
  static vector<int64_t> multiplyPolynomialByNumber(const vector<int64_t>& H, int64_t h, int64_t mod);

  // Function to compute Lagrange basis polynomial L_i(x)
  static vector<int64_t> LagrangePolynomial(int64_t i, const vector<int64_t>& x_values, int64_t mod);

  // Function to compute Lagrange polynomial(x, y)
  static vector<int64_t> setupLagrangePolynomial(const vector<int64_t> x_values, const vector<int64_t> y_values, int64_t mod, const String& name);

  // Function to parse the polynomial string and evaluate it
  static int64_t evaluatePolynomial(const vector<int64_t>& polynomial, int64_t x, int64_t mod);

  // Function to compute the sum of polynomial evaluations at multiple points
  static int64_t sumOfEvaluations(const vector<int64_t>& poly, const vector<int64_t>& points, int64_t mod);

  // Function to create a polynomial for (x - root)
  static vector<int64_t> createLinearPolynomial(int64_t root);

  // Function to calculate Polynomial r(α,x) = (alpha^n - x^n) / (alpha - x)
  static vector<int64_t> calculatePolynomial_r_alpha_x(int64_t alpha, int64_t n, int64_t mod);

  // Function to calculate Polynomial r(α,x) = (alpha^n - x^n) / (alpha - x)
  static int64_t calculatePolynomial_r_alpha_k(int64_t alpha, int64_t k, int64_t n, int64_t mod);

  // Function to expand polynomials given the roots
  static vector<int64_t> expandPolynomials(const vector<int64_t>& roots, int64_t mod);

  // Function to print polynomial in serial
  static void printPolynomial(const vector<int64_t>& coefficients, const String& name);

  // Utility functions for trimming
  static String trim(const String& str);

  // Utility functions for removing commas
  static String removeCommas(const String& str);

  // Utility functions to print a Matrix
  static void printMatrix(vector<vector<int64_t>>& matrix, const String& name);

  // Function to get the row indices of non-zero entries in matrix
  static vector<vector<int64_t>> getNonZeroRows(const vector<vector<int64_t>>& matrix);

  // Function to get the col indices of non-zero entries in matrix
  static vector<vector<int64_t>> getNonZeroCols(const vector<vector<int64_t>>& matrix);

  // Function to create the mapping
  static vector<vector<int64_t>> createMapping(const vector<int64_t>& K, const vector<int64_t>& H, const vector<vector<int64_t>>& nonZero);

  // Function to print the mapping
  static void printMapping(vector<vector<int64_t>>& row, const String& name);

  // Function to create the val mapping
  static vector<vector<int64_t>> valMapping(const vector<int64_t>& K, const vector<int64_t>& H, vector<vector<int64_t>>& nonZeroRows, vector<vector<int64_t>>& nonZeroCols, int64_t mod);

  // Function to calculate log in mod
  static int64_t log_mod(int64_t a, int64_t b, int64_t mod);

  // Function to calculate e_func in mod
  static int64_t e_func(int64_t a, int64_t b, int64_t g, int64_t mod);

  // Function to calculate KZG in mod
  static int64_t KZG_Commitment(vector<int64_t> a, vector<int64_t> b, int64_t mod);
};

#endif  // POLYNOMIAL_H
