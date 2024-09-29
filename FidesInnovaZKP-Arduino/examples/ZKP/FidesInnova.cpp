#include "FidesInnova.h"


// Utility functions for trimming and removing commas
String trim(const String& str) {
  int start = 0;
  while (start < str.length() && isspace(str[start])) start++;
  int end = str.length() - 1;
  while (end >= 0 && isspace(str[end])) end--;
  return str.substring(start, end + 1);
}

String removeCommas(const String& str) {
  String result = str;
  result.replace(",", "");
  return result;
}

void printMatrix(vector<vector<int64_t>>& matrix, const String& name) {
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
vector<vector<int64_t>> getNonZeroRows(const vector<vector<int64_t>>& matrix) {
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
vector<vector<int64_t>> getNonZeroCols(const vector<vector<int64_t>>& matrix) {
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
vector<vector<int64_t>> createMapping(const vector<int64_t>& K, const vector<int64_t>& H, const vector<vector<int64_t>>& nonZero) {
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
void printMapping(vector<vector<int64_t>>& row, const String& name) {
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
vector<vector<int64_t>> valMapping(const vector<int64_t>& K, const vector<int64_t>& H, vector<vector<int64_t>>& nonZeroRows, vector<vector<int64_t>>& nonZeroCols, int64_t mod) {
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














// void FidesInnova::Setup(int64_t g, int64_t d, int64_t l, int64_t mod) {
//   vector<int64_t> pp;
//   int64_t modMinusOne = mod - 1;
//   int64_t current_exponent = 1;
//   int64_t newPower = d % modMinusOne;
//   int64_t res = 0;

//   res = Polynomial::power(g, newPower, mod);
//   pp.push_back(res);
//   for (int64_t i = 1; i < l; ++i) {
//     res = Polynomial::power(res, newPower, mod);
//     pp.push_back(res);
//   }
//   Serial.print("pp = {");
//   for (int64_t i = 0; i < pp.size(); i++) {
//     Serial.print(String(pp[i]));
//     if (pp.size() - i > 1) {
//       Serial.print(", ");
//     }
//   }
//   Serial.println("}");
// }


void FidesInnova::Commitment(String path, int64_t g, int64_t b, int64_t mod) {
  const char* defaultInstructions =
    "li R1, 4\n"
    "mul R1, R1, 5\n"
    "add R1, R1, 11\n"
    "mul R1, R1, 26\n";

  vector<String> instructions;
  // Convert the raw instructions into separate lines
  String instruction(defaultInstructions);
  int startIndex = 0;
  int endIndex;
  // Split the string based on newline character
  while ((endIndex = instruction.indexOf('\n', startIndex)) != -1) {
    instructions.push_back(instruction.substring(startIndex, endIndex));
    startIndex = endIndex + 1;
  }
  // Add the last instruction if there's no newline at the end
  if (startIndex < instruction.length()) {
    instructions.push_back(instruction.substring(startIndex));
  }

  // Number of gates and inputs
  int64_t n_g = 0;
  int64_t n_i = 0;

  // Parse instructions to determine n_g and n_i
  int64_t inputs[32] = { 0 };
  for (const auto& inst : instructions) {
    // Tokenize the instruction manually
    int firstSpace = inst.indexOf(' ');
    String operation = inst.substring(0, firstSpace);   // Extract the operation
    String remainder = inst.substring(firstSpace + 1);  // Extract the remainder (registers and values)

    int secondSpace = remainder.indexOf(' ');
    String rStr = remainder.substring(0, secondSpace);   // Extract first register or value
    String rest = remainder.substring(secondSpace + 1);  // The rest of the instruction

    if (operation == "li") {
      n_i++;
      String xStr = rest;  // Only one value left for li

      int64_t R = rStr.substring(1).toInt();  // Get the register number
      int64_t X = xStr.toInt();               // Get the immediate value

      inputs[R] = X;  // Store the immediate value in the corresponding register
    } else {
      n_g++;

      int thirdSpace = rest.indexOf(' ');
      String xStr = rest.substring(0, thirdSpace);   // Extract the second value (xStr)
      String yStr = rest.substring(thirdSpace + 1);  // Extract the third value (yStr)

      xStr = trim(removeCommas(xStr));
      yStr = trim(removeCommas(yStr));
    }
  }
  Serial.print("Number of immediate instructions (n_i): ");
  Serial.println(n_i);
  Serial.print("Number of general instructions (n_g): ");
  Serial.println(n_g);



  // Matrix order
  int64_t n = n_g + n_i + 1;
  Serial.print("Matrix order: ");
  Serial.println(n);

  // Initialize matrices A, B, C
  vector<vector<int64_t>> A(n, vector<int64_t>(n, 0ll));
  vector<vector<int64_t>> B(n, vector<int64_t>(n, 0ll));
  vector<vector<int64_t>> C(n, vector<int64_t>(n, 0ll));
  // Fill matrices based on the instructions
  int64_t gateIndex = 0;
  int64_t z[n + 1];
  z[0] = 1;

  for (int64_t i = 0; i < n - 1; i++) {
    String inst = instructions[i];
    // Tokenize the instruction
    int firstSpace = inst.indexOf(' ');
    String operation = inst.substring(0, firstSpace);
    String remainder = inst.substring(firstSpace + 1);

    int secondSpace = remainder.indexOf(' ');
    String rStr = remainder.substring(0, secondSpace);
    String rest = remainder.substring(secondSpace + 1);

    int64_t R = rStr.substring(1).toInt();  // Extract register number

    String xStr, yStr;
    if (operation == "li") {
      xStr = rest;
    } else {
      int thirdSpace = rest.indexOf(' ');
      xStr = rest.substring(0, thirdSpace);
      yStr = rest.substring(thirdSpace + 1);
    }

    // Remove commas and trim spaces
    xStr = trim(removeCommas(xStr));
    yStr = trim(removeCommas(yStr));

    int64_t X, Y;

    if (operation == "li") {
      X = xStr.toInt();
      z[i + 1] = X % mod;
    } else {
      gateIndex++;
      int64_t newI = gateIndex + n_i;
      C[newI][newI] = 1;

      if (operation == "add") {
        A[n_i + newI - 1][0] = 1;

        if (isDigit(xStr[0])) {
          X = xStr.toInt();
          z[i + 1] = (z[i] + X) % mod;
          B[n_i + newI - 1][0] = X;
          B[n_i + newI - 1][newI - 1] = 1;
        } else if (isDigit(yStr[0])) {
          Y = yStr.toInt();
          z[i + 1] = (z[i] + Y) % mod;
          B[n_i + newI - 1][0] = Y;
          B[n_i + newI - 1][newI - 1] = 1;
        }
      } else if (operation == "mul") {
        if (isDigit(xStr[0])) {
          X = xStr.toInt();
          z[i + 1] = (z[i] * X) % mod;
          A[n_i + newI - 1][newI - 1] = X;
          B[n_i + newI - 1][newI - 1] = 1;
        } else if (isDigit(yStr[0])) {
          Y = yStr.toInt();
          z[i + 1] = (z[i] * Y) % mod;
          A[n_i + newI - 1][newI - 1] = 1;
          B[n_i + newI - 1][0] = Y;
        }
      }
    }
  }
  Serial.print("z = [");
  for (int i = 0; i < n; i++) {
    Serial.print(z[i]);
    if (n - i > 1) {
      Serial.print(", ");
    }
  }
  Serial.println("]");

  printMatrix(A, "A");
  printMatrix(B, "B");
  printMatrix(C, "C");



  vector<int64_t> H;
  int64_t w, g_n;

  H.push_back(1);
  g_n = ((mod - 1) / n) % mod;
  w = Polynomial::power(g, g_n, mod);
  for (int64_t i = 1; i < n; i++) {
    H.push_back(Polynomial::power(w, i, mod));
  }
  Serial.print("H[n]: ");
  for (int64_t i = 0; i < n; i++) {
    Serial.print(H[i]);
    Serial.print(" ");
  }
  Serial.println("");

  int64_t y, m, t, g_m;

  t = n_i + 1;
  m = (((Polynomial::power(n, 2, mod) - n) / 2) - ((Polynomial::power(t, 2, mod) - t) / 2)) % mod;

  vector<int64_t> K;
  K.push_back(1);
  g_m = ((mod - 1) * Polynomial::modInverse(m, mod)) % mod;
  y = Polynomial::power(g, g_m, mod);
  for (int64_t i = 1; i < m; i++) {
    K.push_back(Polynomial::power(y, i, mod));
  }
  Serial.print("K[m]: ");
  for (int64_t i = 0; i < m; i++) {
    Serial.print(K[i]);
    Serial.print(" ");
  }
  Serial.println("");


  vector<vector<int64_t>> Az(n, vector<int64_t>(1, 0));
  vector<vector<int64_t>> Bz(n, vector<int64_t>(1, 0));
  vector<vector<int64_t>> Cz(n, vector<int64_t>(1, 0));
  // Matrix multiplication with modulo
  for (int64_t i = 0; i < n; i++) {
    for (int64_t j = 0; j < 1; j++) {
      for (int64_t k = 0; k < n; k++) {
        Az[i][j] = (Az[i][j] + (A[i][k] * z[k]) % mod) % mod;
        Bz[i][j] = (Bz[i][j] + (B[i][k] * z[k]) % mod) % mod;
        Cz[i][j] = (Cz[i][j] + (C[i][k] * z[k]) % mod) % mod;
      }
    }
  }
  Serial.print("Matrice Az under modulo ");
  Serial.print(mod);
  Serial.print(" is: ");
  for (int64_t i = 0; i < n; i++) {
    Serial.print(Az[i][0]);
    Serial.print(" ");
  }
  Serial.println("");

  Serial.print("Matrice Bz under modulo ");
  Serial.print(mod);
  Serial.print(" is: ");
  for (int64_t i = 0; i < n; i++) {
    Serial.print(Bz[i][0]);
    Serial.print(" ");
  }
  Serial.println("");

  Serial.print("Matrice Cz under modulo ");
  Serial.print(mod);
  Serial.print(" is: ");
  for (int64_t i = 0; i < n; i++) {
    Serial.print(Cz[i][0]);
    Serial.print(" ");
  }
  Serial.println("");



  vector<vector<int64_t>> zA(2);
  Serial.println("zA(x):");
  for (int64_t i = 0; i < n + b; i++) {
    if (i < n) {
      zA[0].push_back(H[i]);
      zA[1].push_back(Az[i][0]);
    } else {
      zA[0].push_back(Polynomial::generateRandomNumber(H, mod - n));
      zA[1].push_back(Polynomial::generateRandomNumber(H, mod - n));
    }
    Serial.print("zA(");
    Serial.print(zA[0][i]);
    Serial.print(")= ");
    Serial.println(zA[1][i]);
  }
  // zA[0].push_back(30);
  // zA[1].push_back(48);
  // zA[0].push_back(81);
  // zA[1].push_back(81);
  vector<int64_t> z_hatA = Polynomial::setupLagrangePolynomial(zA[0], zA[1], mod, "z_hatA(x)");


  vector<vector<int64_t>> zB(2);
  Serial.println("zB(x):");
  for (int64_t i = 0; i < n + b; i++) {
    if (i < n) {
      zB[0].push_back(H[i]);
      zB[1].push_back(Bz[i][0]);
    } else {
      zB[0].push_back(zA[0][i]);
      zB[1].push_back(Polynomial::generateRandomNumber(H, mod - n));
    }
    Serial.print("zB(");
    Serial.print(zB[0][i]);
    Serial.print(")= ");
    Serial.println(zB[1][i]);
  }
  // zB[0].push_back(30);
  // zB[1].push_back(142);
  // zB[0].push_back(81);
  // zB[1].push_back(156);
  vector<int64_t> z_hatB = Polynomial::setupLagrangePolynomial(zB[0], zB[1], mod, "z_hatB(x)");

  vector<vector<int64_t>> zC(2);
  Serial.println("zC(x):");
  for (int64_t i = 0; i < n + b; i++) {
    if (i < n) {
      zC[0].push_back(H[i]);
      zC[1].push_back(Cz[i][0]);
    } else {
      zC[0].push_back(zA[0][i]);
      zC[1].push_back(Polynomial::generateRandomNumber(H, mod - n));
    }
    Serial.print("zC(");
    Serial.print(zC[0][i]);
    Serial.print(")= ");
    Serial.println(zC[1][i]);
  }
  // zC[0].push_back(30);
  // zC[1].push_back(12);
  // zC[0].push_back(81);
  // zC[1].push_back(129);
  vector<int64_t> z_hatC = Polynomial::setupLagrangePolynomial(zC[0], zC[1], mod, "z_hatC(x)");


  vector<int64_t> zero_to_t_for_H;
  vector<int64_t> t_to_n_for_H;
  vector<int64_t> zero_to_t_for_z;
  vector<int64_t> t_to_n_for_z;
  for (int64_t i = 0; i < t; i++) {
    zero_to_t_for_H.push_back(H[i]);
    zero_to_t_for_z.push_back(z[i]);
  }
  vector<int64_t> polyX_HAT_H = Polynomial::setupLagrangePolynomial(zero_to_t_for_H, zero_to_t_for_z, mod, "x_hat(h)");

  for (int64_t i = 0; i < n - t; i++) {
    t_to_n_for_H.push_back(H[i + t]);
    t_to_n_for_z.push_back(z[i + t]);
  }
  Serial.println("w_bar(h):");
  vector<int64_t> w_bar(n - t + b);
  vector<int64_t> w_bar_numerator(n - t, 1);
  vector<int64_t> w_bar_denominator(n - t, 1);
  for (int64_t i = 0; i < n - t; i++) {
    w_bar_numerator[i] = (t_to_n_for_z[i] - (Polynomial::evaluatePolynomial(polyX_HAT_H, t_to_n_for_H[i], mod)));

    for (int64_t j = 0; j < zero_to_t_for_H.size(); j++) {
      w_bar_denominator[i] *= (t_to_n_for_H[i] - zero_to_t_for_H[j]);
      // Apply modulus to keep the number within the bounds
      w_bar_denominator[i] %= mod;

      // // If result becomes negative, convert it to a positive equivalent under modulo
      if (w_bar_denominator[i] < 0) {
        w_bar_denominator[i] += mod;
      }
    }
    w_bar_denominator[i] = Polynomial::modInverse(w_bar_denominator[i], mod);
    w_bar[i] = (w_bar_numerator[i] * w_bar_denominator[i]) % mod;

    if (w_bar[i] < 0) {
      w_bar[i] += mod;
    }
    Serial.print("w_bar(");
    Serial.print(t_to_n_for_H[i]);
    Serial.print(")= ");
    Serial.println(w_bar[i]);
  }

  Serial.println("w_hat(x):");
  vector<vector<int64_t>> w_hat(2);
  for (int64_t i = 0; i < n - t; i++) {
    w_hat[0].push_back(t_to_n_for_H[i]);
    w_hat[1].push_back(w_bar[i]);
    Serial.print("w_hat(");
    Serial.print(w_hat[0][i]);
    Serial.print(")= ");
    Serial.println(w_hat[1][i]);
  }
  // // w_hat[0].push_back(150);
  // w_hat[1].push_back(42);
  // // w_hat[0].push_back(80);
  // w_hat[1].push_back(180);
  for (int64_t i = n; i < n + b; i++) {
    w_hat[0].push_back(zA[0][i]);
    w_hat[1].push_back(Polynomial::generateRandomNumber(H, mod));
    Serial.print("w_hat(");
    Serial.print(w_hat[0][i - b]);
    Serial.print(")= ");
    Serial.println(w_hat[1][i - b]);
  }
  vector<int64_t> w_hat_x = Polynomial::setupLagrangePolynomial(w_hat[0], w_hat[1], mod, "w_hat(x)");

  vector<int64_t> productAB = Polynomial::multiplyPolynomials(z_hatA, z_hatB, mod);
  vector<int64_t> zAzB_zC = Polynomial::subtractPolynomials(productAB, z_hatC, mod);
  Polynomial::printPolynomial(zAzB_zC, "zA(x)zB(x)-zC(x)");


  vector<int64_t> vH_x(n + 1, 0);
  vH_x[0] = (-1) % mod;
  if (vH_x[0] < 0) {
    vH_x[0] += mod;
  }
  vH_x[n] = 1;
  Polynomial::printPolynomial(vH_x, "vH(x)");


  vector<int64_t> vK_x = Polynomial::createLinearPolynomial(K[0]);
  // Multiply (x - K) for all other K
  for (size_t i = 1; i < K.size(); i++) {
    vector<int64_t> nextPoly = Polynomial::createLinearPolynomial(K[i]);
    vK_x = Polynomial::multiplyPolynomials(vK_x, nextPoly, mod);
  }
  Polynomial::printPolynomial(vK_x, "vK(x)");

  // Dividing the product of zAzB_zC by vH_x
  vector<int64_t> h_0_x = Polynomial::dividePolynomials(zAzB_zC, vH_x, mod)[0];
  Polynomial::printPolynomial(h_0_x, "h0(x)");

  // vector<int64_t> s_x = generateRandomPolynomial(n, (2*n)+b-1, mod);
  vector<int64_t> s_x = { 115, 3, 0, 0, 20, 1, 0, 17, 101, 0, 5 };
  Polynomial::printPolynomial(s_x, "s(x)");

  int64_t sigma1 = Polynomial::sumOfEvaluations(s_x, H, mod);
  Serial.print("sigma1 = ");
  Serial.println(sigma1);

  int64_t alpha = 10;
  int64_t beta1 = 22;  // int64_t beta1 = hashAndExtractLower4Bytes(s_x[8], mod);
  int64_t beta2 = 80;  // int64_t beta2 = hashAndExtractLower4Bytes(s_x[9], mod);
  int64_t etaA = 2;
  int64_t etaB = 30;
  int64_t etaC = 100;

  Serial.print("alpha = ");
  Serial.println(alpha);
  Serial.print("beta1 = ");
  Serial.println(beta1);
  Serial.print("beta2 = ");
  Serial.println(beta2);
  Serial.print("etaA = ");
  Serial.println(etaA);
  Serial.print("etaB = ");
  Serial.println(etaB);
  Serial.print("etaC = ");
  Serial.println(etaC);


  vector<int64_t> etaA_z_hatA_x = Polynomial::multiplyPolynomialByNumber(z_hatA, etaA, mod);
  vector<int64_t> etaB_z_hatB_x = Polynomial::multiplyPolynomialByNumber(z_hatB, etaB, mod);
  vector<int64_t> etaC_z_hatC_x = Polynomial::multiplyPolynomialByNumber(z_hatC, etaC, mod);
  Polynomial::printPolynomial(etaA_z_hatA_x, "etaA_z_hatA(x)");
  Polynomial::printPolynomial(etaB_z_hatB_x, "etaB_z_hatB(x)");
  Polynomial::printPolynomial(etaC_z_hatC_x, "etaC_z_hatC(x)");

  vector<int64_t> Sum_M_eta_M_z_hat_M_x = Polynomial::addPolynomials(Polynomial::addPolynomials(etaA_z_hatA_x, etaB_z_hatB_x, mod), etaC_z_hatC_x, mod);
  Polynomial::printPolynomial(Sum_M_eta_M_z_hat_M_x, "Sum_M_z_hatM(x)");

  vector<int64_t> r_alpha_x = Polynomial::calculatePolynomial_r_alpha_x(alpha, n, mod);
  Polynomial::printPolynomial(r_alpha_x, "r(alpha, x)");

  vector<int64_t> r_Sum_x = Polynomial::multiplyPolynomials(r_alpha_x, Sum_M_eta_M_z_hat_M_x, mod);
  Polynomial::printPolynomial(r_Sum_x, "r(alpha, x)Sum_M_z_hatM(x)");

  vector<int64_t> v_H = Polynomial::expandPolynomials(zero_to_t_for_H, mod);
  Polynomial::printPolynomial(v_H, "v_H");
  vector<int64_t> z_hat_x = Polynomial::addPolynomials(Polynomial::multiplyPolynomials(w_hat_x, v_H, mod), polyX_HAT_H, mod);
  Polynomial::printPolynomial(z_hat_x, "z_hat(x)");



  vector<vector<int64_t>> nonZeroRowsA = getNonZeroRows(A);
  vector<vector<int64_t>> rowA = createMapping(K, H, nonZeroRowsA);
  // rowA[1].push_back(1);
  // rowA[1].push_back(135);
  // rowA[1].push_back(125);
  // rowA[1].push_back(59);
  // rowA[1].push_back(42);
  // rowA[1].push_back(1);
  printMapping(rowA, "row_A");
  vector<vector<int64_t>> nonZeroColsA = getNonZeroCols(A);
  vector<vector<int64_t>> colA = createMapping(K, H, nonZeroColsA);
  // colA[1].push_back(42);
  // colA[1].push_back(1);
  // colA[1].push_back(135);
  // colA[1].push_back(125);
  // colA[1].push_back(59);
  // colA[1].push_back(42);
  printMapping(colA, "col_A");
  vector<vector<int64_t>> valA = valMapping(K, H, nonZeroRowsA, nonZeroColsA, mod);
  printMapping(valA, "val_A");

  vector<int64_t> A_hat(2);
  for (int64_t i = 0; i < nonZeroRowsA[0].size(); i++) {
    vector<int64_t> buff_n = Polynomial::calculatePolynomial_r_alpha_x(rowA[1][i], H.size(), mod);
    int64_t eval = Polynomial::evaluatePolynomial(buff_n, rowA[1][i], mod);
    vector<int64_t> buff = Polynomial::calculatePolynomial_r_alpha_x(colA[1][i], H.size(), mod);
    eval = (eval * valA[1][i]) % mod;
    if (eval < 0) {
      eval += mod;
    }
    buff = Polynomial::multiplyPolynomialByNumber(buff, eval, mod);
    buff = Polynomial::multiplyPolynomialByNumber(buff, Polynomial::calculatePolynomial_r_alpha_k(alpha, rowA[1][i], H.size(), mod), mod);
    if (i > 0) {
      A_hat = Polynomial::addPolynomials(A_hat, buff, mod);
    } else {
      A_hat = buff;
    }
  }
  Polynomial::printPolynomial(A_hat, "A_hat(x)");
  Serial.println("");



  vector<vector<int64_t>> nonZeroRowsB = getNonZeroRows(B);
  vector<vector<int64_t>> rowB = createMapping(K, H, nonZeroRowsB);
  // rowB[1].push_back(59);
  // rowB[1].push_back(1);
  // rowB[1].push_back(42);
  // rowB[1].push_back(135);
  // rowB[1].push_back(59);
  printMapping(rowB, "row_B");
  vector<vector<int64_t>> nonZeroColsB = getNonZeroCols(B);
  vector<vector<int64_t>> colB = createMapping(K, H, nonZeroColsB);
  // colB[1].push_back(59);
  // colB[1].push_back(42);
  // colB[1].push_back(125);
  // colB[1].push_back(1);
  // colB[1].push_back(135);
  printMapping(colB, "col_B");
  vector<vector<int64_t>> valB = valMapping(K, H, nonZeroRowsB, nonZeroColsB, mod);
  printMapping(valB, "val_B");

  vector<int64_t> B_hat(2);
  for (int64_t i = 0; i < nonZeroRowsB[0].size(); i++) {
    vector<int64_t> buff_n = Polynomial::calculatePolynomial_r_alpha_x(rowB[1][i], H.size(), mod);
    int64_t eval = Polynomial::evaluatePolynomial(buff_n, rowB[1][i], mod);
    vector<int64_t> buff = Polynomial::calculatePolynomial_r_alpha_x(colB[1][i], H.size(), mod);
    eval = (eval * valB[1][i]) % mod;
    if (eval < 0) {
      eval += mod;
    }
    buff = Polynomial::multiplyPolynomialByNumber(buff, eval, mod);
    buff = Polynomial::multiplyPolynomialByNumber(buff, Polynomial::calculatePolynomial_r_alpha_k(alpha, rowB[1][i], H.size(), mod), mod);
    if (i > 0) {
      B_hat = Polynomial::addPolynomials(B_hat, buff, mod);
    } else {
      B_hat = buff;
    }
  }
  Polynomial::printPolynomial(B_hat, "B_hat(x)");
  Serial.println("");


  vector<vector<int64_t>> nonZeroRowsC = getNonZeroRows(C);
  vector<vector<int64_t>> rowC = createMapping(K, H, nonZeroRowsC);
  // rowC[1].push_back(1);
  // rowC[1].push_back(59);
  // rowC[1].push_back(125);
  // rowC[1].push_back(1);
  // rowC[1].push_back(135);
  // rowC[1].push_back(42);
  printMapping(rowC, "row_C");
  vector<vector<int64_t>> nonZeroColsC = getNonZeroCols(C);
  vector<vector<int64_t>> colC = createMapping(K, H, nonZeroColsC);
  // colC[1].push_back(125);
  // colC[1].push_back(59);
  // colC[1].push_back(1);
  // colC[1].push_back(1);
  // colC[1].push_back(42);
  // colC[1].push_back(59);
  printMapping(colC, "col_C");
  vector<vector<int64_t>> valC = valMapping(K, H, nonZeroRowsC, nonZeroColsC, mod);
  printMapping(valC, "val_C");

  vector<int64_t> C_hat(2);
  for (int64_t i = 0; i < nonZeroRowsC[0].size(); i++) {
    vector<int64_t> buff_n = Polynomial::calculatePolynomial_r_alpha_x(rowC[1][i], H.size(), mod);
    int64_t eval = Polynomial::evaluatePolynomial(buff_n, rowC[1][i], mod);
    vector<int64_t> buff = Polynomial::calculatePolynomial_r_alpha_x(colC[1][i], H.size(), mod);
    eval = (eval * valC[1][i]) % mod;
    if (eval < 0) {
      eval += mod;
    }
    buff = Polynomial::multiplyPolynomialByNumber(buff, eval, mod);
    buff = Polynomial::multiplyPolynomialByNumber(buff, Polynomial::calculatePolynomial_r_alpha_k(alpha, rowC[1][i], H.size(), mod), mod);
    if (i > 0) {
      C_hat = Polynomial::addPolynomials(C_hat, buff, mod);
    } else {
      C_hat = buff;
    }
  }
  Polynomial::printPolynomial(C_hat, "C_hat(x)");
  Serial.println("");


  vector<int64_t> eta_A_hat = Polynomial::multiplyPolynomialByNumber(A_hat, etaA, mod);
  vector<int64_t> eta_B_hat = Polynomial::multiplyPolynomialByNumber(B_hat, etaB, mod);
  vector<int64_t> eta_C_hat = Polynomial::multiplyPolynomialByNumber(C_hat, etaC, mod);
  Polynomial::printPolynomial(eta_A_hat, "eta_A_hat: ");
  Polynomial::printPolynomial(eta_B_hat, "eta_B_hat: ");
  Polynomial::printPolynomial(eta_C_hat, "eta_C_hat: ");

  vector<int64_t> Sum_M_eta_M_r_M_alpha_x = Polynomial::addPolynomials(Polynomial::addPolynomials(eta_A_hat, eta_B_hat, mod), eta_C_hat, mod);
  Polynomial::printPolynomial(Sum_M_eta_M_r_M_alpha_x, "Sum_M_eta_M_r_M(alpha ,x)");

  vector<int64_t> Sum_M_eta_M_r_M_alpha_x_z_hat_x = Polynomial::multiplyPolynomials(Sum_M_eta_M_r_M_alpha_x, z_hat_x, mod);
  Polynomial::printPolynomial(Sum_M_eta_M_r_M_alpha_x_z_hat_x, "Sum_M_eta_M_r_M(alpha ,x)z-hat(x)");

  vector<int64_t> Sum_check_protocol = Polynomial::addPolynomials(s_x, (Polynomial::subtractPolynomials(r_Sum_x, Sum_M_eta_M_r_M_alpha_x_z_hat_x, mod)), mod);
  Polynomial::printPolynomial(Sum_check_protocol, "Sum_check_protocol");

  vector<int64_t> h_1_x = Polynomial::dividePolynomials(Sum_check_protocol, vH_x, mod)[0];
  Polynomial::printPolynomial(h_1_x, "h1(x)");

  vector<int64_t> g_1_x = Polynomial::dividePolynomials(Sum_check_protocol, vH_x, mod)[1];
  g_1_x.erase(g_1_x.begin());
  Polynomial::printPolynomial(g_1_x, "g1(x)");

  int64_t sigma2 = ((etaA * Polynomial::evaluatePolynomial(A_hat, beta1, mod)) % mod + (etaB * Polynomial::evaluatePolynomial(B_hat, beta1, mod)) % mod + (etaC * Polynomial::evaluatePolynomial(C_hat, beta1, mod)) % mod) % mod;
  Serial.print("sigma2 = ");
  Serial.println(sigma2);


  vector<int64_t> A_hat_M_hat(H.size(), 0);
  vector<int64_t> B_hat_M_hat(H.size(), 0);
  vector<int64_t> C_hat_M_hat(H.size(), 0);
  for (int64_t i = 0; i < nonZeroRowsA[0].size(); i++) {
    vector<int64_t> buff_nA = Polynomial::calculatePolynomial_r_alpha_x(colA[1][i], H.size(), mod);
    int64_t evalA = Polynomial::evaluatePolynomial(buff_nA, beta1, mod);
    vector<int64_t> buffA = Polynomial::calculatePolynomial_r_alpha_x(rowA[1][i], H.size(), mod);
    evalA = (evalA * valA[1][i]) % mod;
    if (evalA < 0) {
      evalA += mod;
    }
    buffA = Polynomial::multiplyPolynomialByNumber(buffA, evalA, mod);

    if (i > 0) {
      A_hat_M_hat = Polynomial::addPolynomials(A_hat_M_hat, buffA, mod);
    } else {
      A_hat_M_hat = buffA;
    }
  }
  for (int64_t i = 0; i < nonZeroRowsB[0].size(); i++) {
    vector<int64_t> buff_nB = Polynomial::calculatePolynomial_r_alpha_x(colB[1][i], H.size(), mod);
    int64_t evalB = Polynomial::evaluatePolynomial(buff_nB, beta1, mod);
    vector<int64_t> buffB = Polynomial::calculatePolynomial_r_alpha_x(rowB[1][i], H.size(), mod);
    evalB = (evalB * valB[1][i]) % mod;
    if (evalB < 0) {
      evalB += mod;
    }
    buffB = Polynomial::multiplyPolynomialByNumber(buffB, evalB, mod);

    if (i > 0) {
      B_hat_M_hat = Polynomial::addPolynomials(B_hat_M_hat, buffB, mod);
    } else {
      B_hat_M_hat = buffB;
    }
  }
  for (int64_t i = 0; i < nonZeroRowsC[0].size(); i++) {
    vector<int64_t> buff_nC = Polynomial::calculatePolynomial_r_alpha_x(colC[1][i], H.size(), mod);
    int64_t evalC = Polynomial::evaluatePolynomial(buff_nC, beta1, mod);
    vector<int64_t> buffC = Polynomial::calculatePolynomial_r_alpha_x(rowC[1][i], H.size(), mod);
    evalC = (evalC * valC[1][i]) % mod;
    if (evalC < 0) {
      evalC += mod;
    }
    buffC = Polynomial::multiplyPolynomialByNumber(buffC, evalC, mod);

    if (i > 0) {
      C_hat_M_hat = Polynomial::addPolynomials(C_hat_M_hat, buffC, mod);
    } else {
      C_hat_M_hat = buffC;
    }
  }
  Polynomial::printPolynomial(A_hat_M_hat, "A_hat_M_hat");
  Polynomial::printPolynomial(B_hat_M_hat, "B_hat_M_hat");
  Polynomial::printPolynomial(C_hat_M_hat, "C_hat_M_hat");

  vector<int64_t> eta_A_hat_M_hat = Polynomial::multiplyPolynomialByNumber(A_hat_M_hat, etaA, mod);
  vector<int64_t> eta_B_hat_M_hat = Polynomial::multiplyPolynomialByNumber(B_hat_M_hat, etaB, mod);
  vector<int64_t> eta_C_hat_M_hat = Polynomial::multiplyPolynomialByNumber(C_hat_M_hat, etaC, mod);
  Polynomial::printPolynomial(eta_A_hat_M_hat, "eta_A_hat_M_hat: ");
  Polynomial::printPolynomial(eta_B_hat_M_hat, "eta_B_hat_M_hat: ");
  Polynomial::printPolynomial(eta_C_hat_M_hat, "eta_C_hat_M_hat: ");

  vector<int64_t> r_Sum_M_eta_M_M_hat_x_beta1 = Polynomial::multiplyPolynomials(Polynomial::addPolynomials(Polynomial::addPolynomials(eta_A_hat_M_hat, eta_B_hat_M_hat, mod), eta_C_hat_M_hat, mod), r_alpha_x, mod);
  Polynomial::printPolynomial(r_Sum_M_eta_M_M_hat_x_beta1, "r_Sum_M_eta_M_M_hat_x_beta1");

  vector<int64_t> h_2_x = Polynomial::dividePolynomials(r_Sum_M_eta_M_M_hat_x_beta1, vH_x, mod)[0];
  Polynomial::printPolynomial(h_2_x, "h2(x)");

  vector<int64_t> g_2_x = Polynomial::dividePolynomials(r_Sum_M_eta_M_M_hat_x_beta1, vH_x, mod)[1];
  g_2_x.erase(g_2_x.begin());
  Polynomial::printPolynomial(g_2_x, "g2(x)");

  vector<int64_t> rowA_x = Polynomial::setupLagrangePolynomial(rowA[0], rowA[1], mod, "rowA(x)");
  vector<int64_t> colA_x = Polynomial::setupLagrangePolynomial(colA[0], colA[1], mod, "colA(x)");
  vector<int64_t> valA_x = Polynomial::setupLagrangePolynomial(valA[0], valA[1], mod, "valA(x)");

  vector<int64_t> rowB_x = Polynomial::setupLagrangePolynomial(rowB[0], rowB[1], mod, "rowB(x)");
  vector<int64_t> colB_x = Polynomial::setupLagrangePolynomial(colB[0], colB[1], mod, "colB(x)");
  vector<int64_t> valB_x = Polynomial::setupLagrangePolynomial(valB[0], valB[1], mod, "valB(x)");

  vector<int64_t> rowC_x = Polynomial::setupLagrangePolynomial(rowC[0], rowC[1], mod, "rowC(x)");
  vector<int64_t> colC_x = Polynomial::setupLagrangePolynomial(colC[0], colC[1], mod, "colC(x)");
  vector<int64_t> valC_x = Polynomial::setupLagrangePolynomial(valC[0], valC[1], mod, "valC(x)");


  int64_t vH_beta1 = Polynomial::evaluatePolynomial(vH_x, beta1, mod);
  Serial.print("vH(beta1) = ");
  Serial.println(vH_beta1);

  int64_t vH_beta2 = Polynomial::evaluatePolynomial(vH_x, beta2, mod);
  Serial.print("vH(beta2) = ");
  Serial.println(vH_beta2);

  vector<int64_t> points_f_3(K.size(), 0);
  int64_t sigma3 = 0;
  for (int64_t i = 0; i < K.size(); i++) {
    int64_t deA = (beta2 - Polynomial::evaluatePolynomial(rowA_x, K[i], mod)) * (beta1 - Polynomial::evaluatePolynomial(colA_x, K[i], mod)) % mod;
    if (deA < 0) deA += mod;

    int64_t deB = (beta2 - Polynomial::evaluatePolynomial(rowB_x, K[i], mod)) * (beta1 - Polynomial::evaluatePolynomial(colB_x, K[i], mod)) % mod;
    if (deB < 0) deB += mod;

    int64_t deC = (beta2 - Polynomial::evaluatePolynomial(rowC_x, K[i], mod)) * (beta1 - Polynomial::evaluatePolynomial(colC_x, K[i], mod)) % mod;
    if (deC < 0) deC += mod;

    int64_t sig3_A = (etaA * (vH_beta2 * vH_beta1 % mod) % mod * Polynomial::evaluatePolynomial(valA_x, K[i], mod) % mod * Polynomial::modInverse(deA, mod)) % mod;
    if (sig3_A < 0) sig3_A += mod;

    int64_t sig3_B = (etaB * (vH_beta2 * vH_beta1 % mod) % mod * Polynomial::evaluatePolynomial(valB_x, K[i], mod) % mod * Polynomial::modInverse(deB, mod)) % mod;
    if (sig3_B < 0) sig3_B += mod;

    int64_t sig3_C = (etaC * (vH_beta2 * vH_beta1 % mod) % mod * Polynomial::evaluatePolynomial(valC_x, K[i], mod) % mod * Polynomial::modInverse(deC, mod)) % mod;
    if (sig3_C < 0) sig3_C += mod;

    points_f_3[i] = (sig3_A + sig3_B + sig3_C) % mod;
    sigma3 += points_f_3[i];
    sigma3 %= mod;
  }
  Serial.print("sigma3 = ");
  Serial.println(sigma3);


  vector<int64_t> poly_beta1 = { beta1 };
  vector<int64_t> poly_beta2 = { beta2 };

  vector<int64_t> poly_pi_a = Polynomial::multiplyPolynomials(Polynomial::subtractPolynomials(rowA_x, poly_beta2, mod), Polynomial::subtractPolynomials(colA_x, poly_beta1, mod), mod);
  vector<int64_t> poly_pi_b = Polynomial::multiplyPolynomials(Polynomial::subtractPolynomials(rowB_x, poly_beta2, mod), Polynomial::subtractPolynomials(colB_x, poly_beta1, mod), mod);
  vector<int64_t> poly_pi_c = Polynomial::multiplyPolynomials(Polynomial::subtractPolynomials(rowC_x, poly_beta2, mod), Polynomial::subtractPolynomials(colC_x, poly_beta1, mod), mod);
  Polynomial::printPolynomial(poly_pi_c, "poly_pi_a");
  Polynomial::printPolynomial(poly_pi_b, "poly_pi_b");
  Polynomial::printPolynomial(poly_pi_c, "poly_pi_c");


  vector<int64_t> poly_etaA_vH_B2_vH_B1 = { (etaA * vH_beta2 * vH_beta1) % mod };
  vector<int64_t> poly_etaB_vH_B2_vH_B1 = { (etaB * vH_beta2 * vH_beta1) % mod };
  vector<int64_t> poly_etaC_vH_B2_vH_B1 = { (etaC * vH_beta2 * vH_beta1) % mod };

  vector<int64_t> poly_sig_a = Polynomial::multiplyPolynomials(poly_etaA_vH_B2_vH_B1, valA_x, mod);
  vector<int64_t> poly_sig_b = Polynomial::multiplyPolynomials(poly_etaB_vH_B2_vH_B1, valB_x, mod);
  vector<int64_t> poly_sig_c = Polynomial::multiplyPolynomials(poly_etaC_vH_B2_vH_B1, valC_x, mod);

  vector<int64_t> a_x = Polynomial::addPolynomials(Polynomial::addPolynomials(Polynomial::multiplyPolynomials(poly_sig_a, Polynomial::multiplyPolynomials(poly_pi_b, poly_pi_c, mod), mod), Polynomial::multiplyPolynomials(poly_sig_b, Polynomial::multiplyPolynomials(poly_pi_a, poly_pi_c, mod), mod), mod), Polynomial::multiplyPolynomials(poly_sig_c, Polynomial::multiplyPolynomials(poly_pi_a, poly_pi_b, mod), mod), mod);
  Polynomial::printPolynomial(a_x, "a(x)");

  vector<int64_t> b_x = Polynomial::multiplyPolynomials(Polynomial::multiplyPolynomials(poly_pi_a, poly_pi_b, mod), poly_pi_c, mod);
  Polynomial::printPolynomial(b_x, "b(x)");


  vector<int64_t> poly_f_3x = Polynomial::setupLagrangePolynomial(K, points_f_3, mod, "poly_f_3(x)");

  vector<int64_t> g_3_x = poly_f_3x;
  g_3_x.erase(g_3_x.begin());
  Polynomial::printPolynomial(g_3_x, "g3(x)");

  vector<int64_t> sigma_3_set_k;
  sigma_3_set_k.push_back((sigma3 * Polynomial::modInverse(K.size(), mod)) % mod);
  Serial.print("sigma_3_set_k = ");
  Serial.println(sigma_3_set_k[0]);

  vector<int64_t> poly_f_3x_new = Polynomial::subtractPolynomials(poly_f_3x, sigma_3_set_k, mod);
  Polynomial::printPolynomial(poly_f_3x_new, "f3(x)new");

  vector<int64_t> h_3_x = Polynomial::dividePolynomials(Polynomial::subtractPolynomials(a_x, Polynomial::multiplyPolynomials(b_x, Polynomial::addPolynomials(poly_f_3x_new, sigma_3_set_k, mod), mod), mod), vK_x, mod)[0];
  Polynomial::printPolynomial(h_3_x, "h3(x)");



  /*
  * Setup
  */
  vector<int64_t> ck;

  // Calculate each expression
  int64_t exp1 = (m) % mod;
  int64_t exp2 = (n-3 + b) % mod;
  int64_t exp3 = (n + b) % mod;
  int64_t exp4 = (n + 2 * b - 1) % mod;
  int64_t exp5 = (2 * n + b - 1) % mod;
  int64_t exp6 = (n + b - 1) % mod;
  int64_t exp7 = (n - 1) % mod;
  int64_t exp8 = (m - 1) % mod;
  int64_t exp9 = (6 * m - 6) % mod;
  
  // Find the maximum value
  int64_t d_AHP = max(exp1, max(exp2, max(exp3, max(exp4, max(exp5, max(exp6, max(exp7, max(exp8, exp9))))))));


  int64_t modMinusOne = mod - 1;
  int64_t current_exponent = 1;
  int64_t newPower = d_AHP % modMinusOne;
  int64_t res = 0;

  // res = Polynomial::power(g, newPower, mod);
  // ck.push_back(res);
  for (int64_t i = 0; i < d_AHP; ++i) {
    res = Polynomial::power(g, Polynomial::power(119, i, mod), mod);
    // res = Polynomial::power(res, newPower, mod);
    ck.push_back(res);
  }
  Serial.print("ck = {");
  for (int64_t i = 0; i < ck.size(); i++) {
    Serial.print(String(ck[i]));
    if (ck.size() - i > 1) {
      Serial.print(", ");
    }
  }
  Serial.println("}");


  int64_t vk = Polynomial::power(g, 119, mod);
  Serial.print("vk = ");
  Serial.println(vk);
  /*
  * Setup
  */




  // }
  // void FidesInnova::Verify() {
  int64_t beta3 = 5;

  bool verify = false;

  int64_t eq11 = (Polynomial::evaluatePolynomial(h_3_x, beta3, mod) * Polynomial::evaluatePolynomial(vK_x, beta3, mod)) % mod;
  int64_t eq12 = (Polynomial::evaluatePolynomial(a_x, beta3, mod) - ((Polynomial::evaluatePolynomial(b_x, beta3, mod) * (beta3 * Polynomial::evaluatePolynomial(g_3_x, beta3, mod) + (sigma3 * Polynomial::modInverse(m, mod)) % mod)))) % mod;
  eq12 %= mod;
  if (eq12 < 0) eq12 += mod;
  Serial.print(eq11);
  Serial.print(" = ");
  Serial.println(eq12);

  int64_t eq21 = (Polynomial::evaluatePolynomial(r_alpha_x, beta2, mod) * sigma3) % mod;
  int64_t eq22 = ((Polynomial::evaluatePolynomial(h_2_x, beta2, mod) * Polynomial::evaluatePolynomial(vH_x, beta2, mod)) % mod + (beta2 * Polynomial::evaluatePolynomial(g_2_x, beta2, mod)) % mod + (sigma2 * Polynomial::modInverse(n, mod)) % mod) % mod;
  eq12 %= mod;
  if (eq12 < 0) eq12 += mod;
  Serial.print(eq21);
  Serial.print(" = ");
  Serial.println(eq22);

  int64_t eq31 = (Polynomial::evaluatePolynomial(s_x, beta1, mod) + Polynomial::evaluatePolynomial(r_alpha_x, beta1, mod) * Polynomial::evaluatePolynomial(Sum_M_eta_M_z_hat_M_x, beta1, mod) - sigma2 * Polynomial::evaluatePolynomial(z_hat_x, beta1, mod));
  int64_t eq32 = (Polynomial::evaluatePolynomial(h_1_x, beta1, mod) * Polynomial::evaluatePolynomial(vH_x, beta1, mod) + beta1 * Polynomial::evaluatePolynomial(g_1_x, beta1, mod) + sigma1 * Polynomial::modInverse(n, mod)) % mod;
  eq31 %= mod;
  if (eq31 < 0) eq31 += mod;
  Serial.print(eq31);
  Serial.print(" = ");
  Serial.println(eq32);

  int64_t eq41 = (Polynomial::evaluatePolynomial(z_hatA, beta1, mod) * Polynomial::evaluatePolynomial(z_hatB, beta1, mod) - Polynomial::evaluatePolynomial(z_hatC, beta1, mod));
  int64_t eq42 = (Polynomial::evaluatePolynomial(h_0_x, beta1, mod) * Polynomial::evaluatePolynomial(vH_x, beta1, mod)) % mod;
  eq41 %= mod;
  if (eq41 < 0) eq41 += mod;

  if (eq11 == eq12 && eq21 == eq22 && eq31 == eq32 && eq41 == eq42) {
    verify = true;
  }
  Serial.print(eq41);
  Serial.print(" = ");
  Serial.println(eq42);


  Serial.println("");
  if (verify) {
    Serial.println("verify!!!!!!!!!!");
  }
  Serial.print("Time: ");
  Serial.print(millis() / 1000);
  Serial.println(" S");
}
