/*
  Read the setup.json and program_param.json from flash of the RISC-V microcontroller

  setup.json
  {
    "Class":  32-bit Integer,
    "ck": 64-bit Integer Array,
    "vk": 64-bit Integer
  }

  
  program_param.json
  {
    "Class":  32-bit Integer,
    "A": 64-bit Integer Array,
    "B": 64-bit Integer Array<64-bit Integer Array>
  }

*/

#include "fidesinnova.h"

void fidesinnova::proofGenerator(String path, int64_t g, int64_t b, int64_t mod) {
  //It's duplicate codes as the commitment
  String setup = readFile("/setup.json");
  DynamicJsonDocument jsonSetup(2048);  // Create a DynamicJsonDocument with a buffer size
  deserializeJson(jsonSetup, setup);
  vector<int64_t> ck;
  JsonArray array = jsonSetup["ck"];
  for (JsonVariant v : array) {
    ck.push_back(v.as<int64_t>());
  }
  int64_t vk = jsonSetup["vk"][0].as<int64_t>();

  const char* defaultInstructions =
    "li R1, 4\n"
    "mul R1, R1, 5\n"
    "addi R1, R1, 11\n"
    "mul R1, R1, 26\n";

  vector<String> instructions;
  // Convert the raw instructions into separate lines
  String instruction(defaultInstructions);
  int64_t startIndex = 0;
  int64_t endIndex;
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
    int64_t firstSpace = inst.indexOf(' ');
    String operation = inst.substring(0, firstSpace);   // Extract the operation
    String remainder = inst.substring(firstSpace + 1);  // Extract the remainder (registers and values)

    int64_t secondSpace = remainder.indexOf(' ');
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

      int64_t thirdSpace = rest.indexOf(' ');
      String xStr = rest.substring(0, thirdSpace);   // Extract the second value (xStr)
      String yStr = rest.substring(thirdSpace + 1);  // Extract the third value (yStr)

      xStr = Polynomial::trim(Polynomial::removeCommas(xStr));
      yStr = Polynomial::trim(Polynomial::removeCommas(yStr));
    }
  }
  Serial.print("Number of immediate instructions (n_i): ");
  Serial.println(n_i);
  Serial.print("Number of general instructions (n_g): ");
  Serial.println(n_g);

  // Matrix order
  int64_t n = n_g + n_i + 1;
  int64_t m, t;
  Serial.print("Matrix order: ");
  Serial.println(n);

  t = n_i + 1;
  m = (((Polynomial::power(n, 2, mod) - n) / 2) - ((Polynomial::power(t, 2, mod) - t) / 2)) % mod;

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
    int64_t firstSpace = inst.indexOf(' ');
    String operation = inst.substring(0, firstSpace);
    String remainder = inst.substring(firstSpace + 1);

    int64_t secondSpace = remainder.indexOf(' ');
    String rStr = remainder.substring(0, secondSpace);
    String rest = remainder.substring(secondSpace + 1);

    int64_t R = rStr.substring(1).toInt();  // Extract register number

    String xStr, yStr;
    if (operation == "li") {
      xStr = rest;
    } else {
      int64_t thirdSpace = rest.indexOf(' ');
      xStr = rest.substring(0, thirdSpace);
      yStr = rest.substring(thirdSpace + 1);
    }

    // Remove commas and trim spaces
    xStr = Polynomial::trim(Polynomial::removeCommas(xStr));
    yStr = Polynomial::trim(Polynomial::removeCommas(yStr));

    int64_t X, Y;

    if (operation == "li") {
      X = xStr.toInt();
      z[i + 1] = X % mod;
    } else {
      gateIndex++;
      int64_t newI = gateIndex + n_i;
      C[newI][newI] = 1;

      if (operation == "addi") {
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
  for (int64_t i = 0; i < n; i++) {
    Serial.print(z[i]);
    if (n - i > 1) {
      Serial.print(", ");
    }
  }
  Serial.println("]");

  Polynomial::printMatrix(A, "A");
  Polynomial::printMatrix(B, "B");
  Polynomial::printMatrix(C, "C");

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

  int64_t y, g_m;

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
  // for (int64_t i = 0; i < n + b; i++) {
  for (int64_t i = 0; i < n; i++) {
    if (i < n) {
      zA[0].push_back(H[i]);
      zA[1].push_back(Az[i][0]);
    } else {
      // zA[0].push_back(Polynomial::generateRandomNumber(H, mod - n));
      // zA[1].push_back(Polynomial::generateRandomNumber(H, mod - n));
    }
    Serial.print("zA(");
    Serial.print(zA[0][i]);
    Serial.print(")= ");
    Serial.println(zA[1][i]);
  }
  zA[0].push_back(150);
  zA[1].push_back(5);
  zA[0].push_back(80);
  zA[1].push_back(47);

  vector<int64_t> z_hatA = Polynomial::setupLagrangePolynomial(zA[0], zA[1], mod, "z_hatA(x)");


  vector<vector<int64_t>> zB(2);
  Serial.println("zB(x):");
  // for (int64_t i = 0; i < n + b; i++) {
  for (int64_t i = 0; i < n; i++) {
    if (i < n) {
      zB[0].push_back(H[i]);
      zB[1].push_back(Bz[i][0]);
    } else {
      // zB[0].push_back(zA[0][i]);
      // zB[1].push_back(Polynomial::generateRandomNumber(H, mod - n));
    }
    Serial.print("zB(");
    Serial.print(zB[0][i]);
    Serial.print(")= ");
    Serial.println(zB[1][i]);
  }
  zB[0].push_back(150);
  zB[1].push_back(15);
  zB[0].push_back(80);
  zB[1].push_back(170);
  vector<int64_t> z_hatB = Polynomial::setupLagrangePolynomial(zB[0], zB[1], mod, "z_hatB(x)");

  vector<vector<int64_t>> zC(2);
  Serial.println("zC(x):");
  // for (int64_t i = 0; i < n + b; i++) {
  for (int64_t i = 0; i < n; i++) {
    if (i < n) {
      zC[0].push_back(H[i]);
      zC[1].push_back(Cz[i][0]);
    } else {
      // zC[0].push_back(zA[0][i]);
      // zC[1].push_back(Polynomial::generateRandomNumber(H, mod - n));
    }
    Serial.print("zC(");
    Serial.print(zC[0][i]);
    Serial.print(")= ");
    Serial.println(zC[1][i]);
  }
  zC[0].push_back(150);
  zC[1].push_back(1);
  zC[0].push_back(80);
  zC[1].push_back(100);
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
  w_hat[0].push_back(150);
  w_hat[1].push_back(42);
  w_hat[0].push_back(80);
  w_hat[1].push_back(180);
  for (int64_t i = n; i < n + b; i++) {
    // w_hat[0].push_back(zA[0][i]);
    // w_hat[1].push_back(Polynomial::generateRandomNumber(H, mod));
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
    vH_x[0] += mod; //Handling negative values properly
  }
  vH_x[n] = 1;
  Polynomial::printPolynomial(vH_x, "KvH(x)");


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

  vector<vector<int64_t>> nonZeroRowsA = Polynomial::getNonZeroRows(A);
  vector<vector<int64_t>> rowA = Polynomial::createMapping(K, H, nonZeroRowsA);
  rowA[1].push_back(1);
  rowA[1].push_back(135);
  rowA[1].push_back(125);
  rowA[1].push_back(59);
  rowA[1].push_back(42);
  rowA[1].push_back(1);
  Polynomial::printMapping(rowA, "row_A");
  vector<vector<int64_t>> nonZeroColsA = Polynomial::getNonZeroCols(A);
  vector<vector<int64_t>> colA = Polynomial::createMapping(K, H, nonZeroColsA);
  colA[1].push_back(42);
  colA[1].push_back(1);
  colA[1].push_back(135);
  colA[1].push_back(125);
  colA[1].push_back(59);
  colA[1].push_back(42);
  Polynomial::printMapping(colA, "col_A");
  vector<vector<int64_t>> valA = Polynomial::valMapping(K, H, nonZeroRowsA, nonZeroColsA, mod);
  Polynomial::printMapping(valA, "val_A");

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

  vector<vector<int64_t>> nonZeroRowsB = Polynomial::getNonZeroRows(B);
  vector<vector<int64_t>> rowB = Polynomial::createMapping(K, H, nonZeroRowsB);
  rowB[1].push_back(59);
  rowB[1].push_back(1);
  rowB[1].push_back(42);
  rowB[1].push_back(135);
  rowB[1].push_back(59);
  Polynomial::printMapping(rowB, "row_B");
  vector<vector<int64_t>> nonZeroColsB = Polynomial::getNonZeroCols(B);
  vector<vector<int64_t>> colB = Polynomial::createMapping(K, H, nonZeroColsB);
  colB[1].push_back(59);
  colB[1].push_back(42);
  colB[1].push_back(125);
  colB[1].push_back(1);
  colB[1].push_back(135);
  Polynomial::printMapping(colB, "col_B");
  vector<vector<int64_t>> valB = Polynomial::valMapping(K, H, nonZeroRowsB, nonZeroColsB, mod);
  Polynomial::printMapping(valB, "val_B");

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


  vector<vector<int64_t>> nonZeroRowsC = Polynomial::getNonZeroRows(C);
  vector<vector<int64_t>> rowC = Polynomial::createMapping(K, H, nonZeroRowsC);
  rowC[1].push_back(1);
  rowC[1].push_back(59);
  rowC[1].push_back(125);
  rowC[1].push_back(1);
  rowC[1].push_back(135);
  rowC[1].push_back(42);
  Polynomial::printMapping(rowC, "row_C");
  vector<vector<int64_t>> nonZeroColsC = Polynomial::getNonZeroCols(C);
  vector<vector<int64_t>> colC = Polynomial::createMapping(K, H, nonZeroColsC);
  colC[1].push_back(125);
  colC[1].push_back(59);
  colC[1].push_back(1);
  colC[1].push_back(1);
  colC[1].push_back(42);
  colC[1].push_back(59);
  Polynomial::printMapping(colC, "col_C");
  vector<vector<int64_t>> valC = Polynomial::valMapping(K, H, nonZeroRowsC, nonZeroColsC, mod);
  Polynomial::printMapping(valC, "val_C");

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

/*
// Check for overflow risk for etaA, etaB and etaC and ensure mod is applied correctly
// vector<int64_t> etaA_z_hatA_x = Polynomial::multiplyPolynomialByNumber(A_hatA, etaA % mod, mod);
// vector<int64_t> etaB_z_hatB_x = Polynomial::multiplyPolynomialByNumber(B_hatB, etaB % mod, mod);
// vector<int64_t> etaC_z_hatC_x = Polynomial::multiplyPolynomialByNumber(C_hatC, etaC % mod, mod);*/

  vector<int64_t> eta_A_hat = Polynomial::multiplyPolynomialByNumber(A_hat, etaA, mod);
  vector<int64_t> eta_B_hat = Polynomial::multiplyPolynomialByNumber(B_hat, etaB, mod);
  vector<int64_t> eta_C_hat = Polynomial::multiplyPolynomialByNumber(C_hat, etaC, mod);
  Polynomial::printPolynomial(eta_A_hat, "eta_A_hat: ");
  Polynomial::printPolynomial(eta_B_hat, "eta_B_hat: ");
  Polynomial::printPolynomial(eta_C_hat, "eta_C_hat: ");

  // Calculate the sum of the three polynomials and print the result
  vector<int64_t> Sum_M_eta_M_r_M_alpha_x = Polynomial::addPolynomials(Polynomial::addPolynomials(eta_A_hat, eta_B_hat, mod), eta_C_hat, mod);
  Polynomial::printPolynomial(Sum_M_eta_M_r_M_alpha_x, "Sum_M_eta_M_r_M(alpha ,x)");

  // Multiply the sum by another polynomial z_hat_x and print the result
  vector<int64_t> Sum_M_eta_M_r_M_alpha_x_z_hat_x = Polynomial::multiplyPolynomials(Sum_M_eta_M_r_M_alpha_x, z_hat_x, mod);
  Polynomial::printPolynomial(Sum_M_eta_M_r_M_alpha_x_z_hat_x, "Sum_M_eta_M_r_M(alpha ,x)z-hat(x)");

  // Calculate the sum for the check protocol, subtracting the modified sum from s_x
  vector<int64_t> Sum_check_protocol = Polynomial::addPolynomials(s_x, (Polynomial::subtractPolynomials(r_Sum_x, Sum_M_eta_M_r_M_alpha_x_z_hat_x, mod)), mod);
  Polynomial::printPolynomial(Sum_check_protocol, "Sum_check_protocol");

  // Divide the sum check protocol by vH_x to get two results: h1(x) and g1(x)
  vector<int64_t> h_1_x = Polynomial::dividePolynomials(Sum_check_protocol, vH_x, mod)[0];
  Polynomial::printPolynomial(h_1_x, "h1(x)");

  // Get the second part of the division result, g1(x), and erase the first element
  vector<int64_t> g_1_x = Polynomial::dividePolynomials(Sum_check_protocol, vH_x, mod)[1];
  g_1_x.erase(g_1_x.begin());
  Polynomial::printPolynomial(g_1_x, "g1(x)");

  // Calculate sigma2 using the evaluations of the polynomials A_hat, B_hat, and C_hat and print the result of sigma2
  int64_t sigma2 = ((etaA * Polynomial::evaluatePolynomial(A_hat, beta1, mod)) % mod + (etaB * Polynomial::evaluatePolynomial(B_hat, beta1, mod)) % mod + (etaC * Polynomial::evaluatePolynomial(C_hat, beta1, mod)) % mod) % mod;
  Serial.print("sigma2 = ");
  Serial.println(sigma2);

  // Initialize vectors for the modified polynomial results with zeros
  vector<int64_t> A_hat_M_hat(H.size(), 0);
  vector<int64_t> B_hat_M_hat(H.size(), 0);
  vector<int64_t> C_hat_M_hat(H.size(), 0);

  // Loop through non-zero rows for matrix A and calculate the modified polynomial A_hat_M_hat
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
  // Print the final modified polynomials for A, B, and C
  Polynomial::printPolynomial(A_hat_M_hat, "A_hat_M_hat");
  Polynomial::printPolynomial(B_hat_M_hat, "B_hat_M_hat");
  Polynomial::printPolynomial(C_hat_M_hat, "C_hat_M_hat");

  // Multiply the modified polynomials by their respective eta values and print
  vector<int64_t> eta_A_hat_M_hat = Polynomial::multiplyPolynomialByNumber(A_hat_M_hat, etaA, mod);
  vector<int64_t> eta_B_hat_M_hat = Polynomial::multiplyPolynomialByNumber(B_hat_M_hat, etaB, mod);
  vector<int64_t> eta_C_hat_M_hat = Polynomial::multiplyPolynomialByNumber(C_hat_M_hat, etaC, mod);
  Polynomial::printPolynomial(eta_A_hat_M_hat, "eta_A_hat_M_hat: ");
  Polynomial::printPolynomial(eta_B_hat_M_hat, "eta_B_hat_M_hat: ");
  Polynomial::printPolynomial(eta_C_hat_M_hat, "eta_C_hat_M_hat: ");

  // Calculate the final result for r_Sum_M_eta_M_M_hat_x_beta1
  vector<int64_t> r_Sum_M_eta_M_M_hat_x_beta1 = Polynomial::multiplyPolynomials(Polynomial::addPolynomials(Polynomial::addPolynomials(eta_A_hat_M_hat, eta_B_hat_M_hat, mod), eta_C_hat_M_hat, mod), r_alpha_x, mod);
  Polynomial::printPolynomial(r_Sum_M_eta_M_M_hat_x_beta1, "r_Sum_M_eta_M_M_hat_x_beta1");

  // Divide the final result by vH_x to get h2(x) and g2(x)
  vector<int64_t> h_2_x = Polynomial::dividePolynomials(r_Sum_M_eta_M_M_hat_x_beta1, vH_x, mod)[0];
  Polynomial::printPolynomial(h_2_x, "h2(x)");

  vector<int64_t> g_2_x = Polynomial::dividePolynomials(r_Sum_M_eta_M_M_hat_x_beta1, vH_x, mod)[1];
  g_2_x.erase(g_2_x.begin());//remove the first item
  Polynomial::printPolynomial(g_2_x, "g2(x)");

  // Generate Lagrange polynomials for row, column, and value vectors for matrices A and B
  vector<int64_t> rowA_x = Polynomial::setupLagrangePolynomial(rowA[0], rowA[1], mod, "rowA(x)");
  vector<int64_t> colA_x = Polynomial::setupLagrangePolynomial(colA[0], colA[1], mod, "colA(x)");
  vector<int64_t> valA_x = Polynomial::setupLagrangePolynomial(valA[0], valA[1], mod, "valA(x)");

  vector<int64_t> rowB_x = Polynomial::setupLagrangePolynomial(rowB[0], rowB[1], mod, "rowB(x)");
  vector<int64_t> colB_x = Polynomial::setupLagrangePolynomial(colB[0], colB[1], mod, "colB(x)");
  vector<int64_t> valB_x = Polynomial::setupLagrangePolynomial(valB[0], valB[1], mod, "valB(x)");

  vector<int64_t> rowC_x = Polynomial::setupLagrangePolynomial(rowC[0], rowC[1], mod, "rowC(x)");
  vector<int64_t> colC_x = Polynomial::setupLagrangePolynomial(colC[0], colC[1], mod, "colC(x)");
  vector<int64_t> valC_x = Polynomial::setupLagrangePolynomial(valC[0], valC[1], mod, "valC(x)");

  // Evaluate polynomial vH at beta1 and beta2
  int64_t vH_beta1 = Polynomial::evaluatePolynomial(vH_x, beta1, mod);
  Serial.print("vH(beta1) = ");
  Serial.println(vH_beta1);

  int64_t vH_beta2 = Polynomial::evaluatePolynomial(vH_x, beta2, mod);
  Serial.print("vH(beta2) = ");
  Serial.println(vH_beta2);

  // Initialize vectors for function points and sigma value
  vector<int64_t> points_f_3(K.size(), 0);
  int64_t sigma3 = 0;

  // Loop over K to compute delta and signature values for A, B, and C
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

  // Create polynomials for beta1 and beta2
  vector<int64_t> poly_beta1 = { beta1 };
  vector<int64_t> poly_beta2 = { beta2 };

  // Compute polynomial products for signatures
  vector<int64_t> poly_pi_a = Polynomial::multiplyPolynomials(Polynomial::subtractPolynomials(rowA_x, poly_beta2, mod), Polynomial::subtractPolynomials(colA_x, poly_beta1, mod), mod);
  vector<int64_t> poly_pi_b = Polynomial::multiplyPolynomials(Polynomial::subtractPolynomials(rowB_x, poly_beta2, mod), Polynomial::subtractPolynomials(colB_x, poly_beta1, mod), mod);
  vector<int64_t> poly_pi_c = Polynomial::multiplyPolynomials(Polynomial::subtractPolynomials(rowC_x, poly_beta2, mod), Polynomial::subtractPolynomials(colC_x, poly_beta1, mod), mod);
  Polynomial::printPolynomial(poly_pi_c, "poly_pi_a");
  Polynomial::printPolynomial(poly_pi_b, "poly_pi_b");
  Polynomial::printPolynomial(poly_pi_c, "poly_pi_c");

  // Compute polynomials for signature multipliers
  vector<int64_t> poly_etaA_vH_B2_vH_B1 = { (etaA * vH_beta2 * vH_beta1) % mod };
  vector<int64_t> poly_etaB_vH_B2_vH_B1 = { (etaB * vH_beta2 * vH_beta1) % mod };
  vector<int64_t> poly_etaC_vH_B2_vH_B1 = { (etaC * vH_beta2 * vH_beta1) % mod };

  // Calculate signatures
  vector<int64_t> poly_sig_a = Polynomial::multiplyPolynomials(poly_etaA_vH_B2_vH_B1, valA_x, mod);
  vector<int64_t> poly_sig_b = Polynomial::multiplyPolynomials(poly_etaB_vH_B2_vH_B1, valB_x, mod);
  vector<int64_t> poly_sig_c = Polynomial::multiplyPolynomials(poly_etaC_vH_B2_vH_B1, valC_x, mod);

  vector<int64_t> a_x = Polynomial::addPolynomials(Polynomial::addPolynomials(Polynomial::multiplyPolynomials(poly_sig_a, Polynomial::multiplyPolynomials(poly_pi_b, poly_pi_c, mod), mod), Polynomial::multiplyPolynomials(poly_sig_b, Polynomial::multiplyPolynomials(poly_pi_a, poly_pi_c, mod), mod), mod), Polynomial::multiplyPolynomials(poly_sig_c, Polynomial::multiplyPolynomials(poly_pi_a, poly_pi_b, mod), mod), mod);
  Polynomial::printPolynomial(a_x, "a(x)");

  vector<int64_t> b_x = Polynomial::multiplyPolynomials(Polynomial::multiplyPolynomials(poly_pi_a, poly_pi_b, mod), poly_pi_c, mod);
  Polynomial::printPolynomial(b_x, "b(x)");

  // Set up polynomial for f_3 using K
  vector<int64_t> poly_f_3x = Polynomial::setupLagrangePolynomial(K, points_f_3, mod, "poly_f_3(x)");

  vector<int64_t> g_3_x = poly_f_3x;
  g_3_x.erase(g_3_x.begin());
  Polynomial::printPolynomial(g_3_x, "g3(x)");

  // Calculate sigma_3_set_k based on sigma3 and K.size()
  vector<int64_t> sigma_3_set_k;
  sigma_3_set_k.push_back((sigma3 * Polynomial::modInverse(K.size(), mod)) % mod);
  Serial.print("sigma_3_set_k = ");
  Serial.println(sigma_3_set_k[0]);

  // Update polynomial f_3 by subtracting sigma_3_set_k
  vector<int64_t> poly_f_3x_new = Polynomial::subtractPolynomials(poly_f_3x, sigma_3_set_k, mod);
  Polynomial::printPolynomial(poly_f_3x_new, "f3(x)new");

  // Calculate polynomial h_3(x) using previous results
  vector<int64_t> h_3_x = Polynomial::dividePolynomials(Polynomial::subtractPolynomials(a_x, Polynomial::multiplyPolynomials(b_x, Polynomial::addPolynomials(poly_f_3x_new, sigma_3_set_k, mod), mod), mod), vK_x, mod)[0];
  Polynomial::printPolynomial(h_3_x, "h3(x)");

  // Define random values based on s_x
  int64_t eta_w_hat = 1;    // a random number  based on s_x
  int64_t eta_z_hatA = 4;   // a random number  based on s_x
  int64_t eta_z_hatB = 10;  // a random number  based on s_x
  int64_t eta_z_hatC = 8;   // a random number  based on s_x
  int64_t eta_h_0_x = 32;   // a random number  based on s_x
  int64_t eta_s_x = 45;     // a random number  based on s_x
  int64_t eta_g_1_x = 92;   // a random number  based on s_x
  int64_t eta_h_1_x = 11;   // a random number  based on s_x
  int64_t eta_g_2_x = 1;    // a random number  based on s_x
  int64_t eta_h_2_x = 5;    // a random number  based on s_x
  int64_t eta_g_3_x = 25;   // a random number  based on s_x
  int64_t eta_h_3_x = 63;   // a random number  based on s_x

  // Initialize the polynomial p(x) by performing several polynomial operations and print
  vector<int64_t> p_x =
    Polynomial::addPolynomials(
      Polynomial::addPolynomials(Polynomial::addPolynomials(Polynomial::multiplyPolynomialByNumber(w_hat_x, eta_w_hat, mod), Polynomial::multiplyPolynomialByNumber(z_hatA, eta_z_hatA, mod), mod),
                                 Polynomial::addPolynomials(Polynomial::multiplyPolynomialByNumber(z_hatB, eta_z_hatB, mod), Polynomial::multiplyPolynomialByNumber(z_hatC, eta_z_hatC, mod), mod), mod),
      Polynomial::addPolynomials(
        Polynomial::addPolynomials(Polynomial::addPolynomials(Polynomial::multiplyPolynomialByNumber(h_0_x, eta_h_0_x, mod), Polynomial::multiplyPolynomialByNumber(s_x, eta_s_x, mod), mod),
                                   Polynomial::addPolynomials(Polynomial::multiplyPolynomialByNumber(g_1_x, eta_g_1_x, mod), Polynomial::multiplyPolynomialByNumber(h_1_x, eta_h_1_x, mod), mod), mod),
        Polynomial::addPolynomials(Polynomial::addPolynomials(Polynomial::multiplyPolynomialByNumber(g_2_x, eta_g_2_x, mod), Polynomial::multiplyPolynomialByNumber(h_2_x, eta_h_2_x, mod), mod),
                                   Polynomial::addPolynomials(Polynomial::multiplyPolynomialByNumber(g_3_x, eta_g_3_x, mod), Polynomial::multiplyPolynomialByNumber(h_3_x, eta_h_3_x, mod), mod), mod),
        mod),
      mod);
  Polynomial::printPolynomial(p_x, "p(x)");

  // Choose a random value for z
  int64_t z_random = 2;
  int64_t y_prime = Polynomial::evaluatePolynomial(p_x, z_random, mod);
  Serial.print("y_prime = ");
  Serial.println(y_prime);

  // p(x) - y'  =>  p(x)  !!!!!!!
  vector<int64_t> q_xBuf;
  q_xBuf.push_back(mod - z_random);
  q_xBuf.push_back(1);
  Polynomial::printPolynomial(q_xBuf, "div = ");

  vector<int64_t> q_x = Polynomial::dividePolynomials(p_x, q_xBuf, mod)[0];
  Polynomial::printPolynomial(q_x, "q(x)");

  // Generate a KZG commitment for q(x) using the provided verification key (ck)
  int64_t p_17_AHP = Polynomial::KZG_Commitment(ck, q_x, mod);
  Serial.print("p_17_AHP = ");
  Serial.println(p_17_AHP);

  // Generate KZG commitments for various polynomials and print
  int64_t Com1_AHP_x = Polynomial::KZG_Commitment(ck, w_hat_x, mod);
  int64_t Com2_AHP_x = Polynomial::KZG_Commitment(ck, z_hatA, mod);
  int64_t Com3_AHP_x = Polynomial::KZG_Commitment(ck, z_hatB, mod);
  int64_t Com4_AHP_x = Polynomial::KZG_Commitment(ck, z_hatC, mod);
  int64_t Com5_AHP_x = Polynomial::KZG_Commitment(ck, h_0_x, mod);
  int64_t Com6_AHP_x = Polynomial::KZG_Commitment(ck, s_x, mod);
  int64_t Com7_AHP_x = Polynomial::KZG_Commitment(ck, g_1_x, mod);
  int64_t Com8_AHP_x = Polynomial::KZG_Commitment(ck, h_1_x, mod);
  int64_t Com9_AHP_x = Polynomial::KZG_Commitment(ck, g_2_x, mod);
  int64_t Com10_AHP_x = Polynomial::KZG_Commitment(ck, h_2_x, mod);
  int64_t Com11_AHP_x = Polynomial::KZG_Commitment(ck, g_3_x, mod);
  int64_t Com12_AHP_x = Polynomial::KZG_Commitment(ck, h_3_x, mod);
  Serial.print("Com1_AHP_x = ");
  Serial.println(Com1_AHP_x);
  Serial.print("Com2_AHP_x = ");
  Serial.println(Com2_AHP_x);
  Serial.print("Com3_AHP_x = ");
  Serial.println(Com3_AHP_x);
  Serial.print("Com4_AHP_x = ");
  Serial.println(Com4_AHP_x);
  Serial.print("Com5_AHP_x = ");
  Serial.println(Com5_AHP_x);
  Serial.print("Com6_AHP_x = ");
  Serial.println(Com6_AHP_x);
  Serial.print("Com7_AHP_x = ");
  Serial.println(Com7_AHP_x);
  Serial.print("Com8_AHP_x = ");
  Serial.println(Com8_AHP_x);
  Serial.print("Com9_AHP_x = ");
  Serial.println(Com9_AHP_x);
  Serial.print("Com10_AHP_x = ");
  Serial.println(Com10_AHP_x);
  Serial.print("Com11_AHP_x = ");
  Serial.println(Com11_AHP_x);
  Serial.print("Com12_AHP_x = ");
  Serial.println(Com12_AHP_x);

  // Generate a KZG commitment for the combined polynomial p(x)
  int64_t ComP_AHP_x = Polynomial::KZG_Commitment(ck, p_x, mod);
  Serial.print("ComP_AHP = ");
  Serial.println(ComP_AHP_x);



  // TODO: Add req params to the proof.json
  /*
    "CommitmentID":  64-bit,
    "Class": 32-bit Integer,
    "IoT_Manufacturer_Name": String,
    "IoT_Device_Name": String,
    "Device_Hardware_Version": float,
    "Firmware_Version": float,
    "Lines": 64-bit Array,
  */



  //create a Json document to store proof-realated data
  DynamicJsonDocument doc(2048);
  JsonArray jsonArray;
  jsonArray = doc.createNestedArray("P1AHP");
  jsonArray.add(sigma1);
  jsonArray = doc.createNestedArray("P2AHP");
  for (int64_t value : w_hat_x) {
    jsonArray.add(value);
  }
  jsonArray = doc.createNestedArray("P3AHP");
  for (int64_t value : z_hatA) {
    jsonArray.add(value);
  }
  jsonArray = doc.createNestedArray("P4AHP");
  for (int64_t value : z_hatB) {
    jsonArray.add(value);
  }
  jsonArray = doc.createNestedArray("P5AHP");
  for (int64_t value : z_hatC) {
    jsonArray.add(value);
  }
  jsonArray = doc.createNestedArray("P6AHP");
  for (int64_t value : h_0_x) {
    jsonArray.add(value);
  }
  jsonArray = doc.createNestedArray("P7AHP");
  for (int64_t value : s_x) {
    jsonArray.add(value);
  }
  jsonArray = doc.createNestedArray("P8AHP");
  for (int64_t value : g_1_x) {
    jsonArray.add(value);
  }
  jsonArray = doc.createNestedArray("P9AHP");
  for (int64_t value : h_1_x) {
    jsonArray.add(value);
  }
  jsonArray = doc.createNestedArray("P10AHP");
  jsonArray.add(sigma2);
  jsonArray = doc.createNestedArray("P11AHP");
  for (int64_t value : g_2_x) {
    jsonArray.add(value);
  }
  jsonArray = doc.createNestedArray("P12AHP");
  for (int64_t value : h_2_x) {
    jsonArray.add(value);
  }
  jsonArray = doc.createNestedArray("P13AHP");
  jsonArray.add(sigma3);
  jsonArray = doc.createNestedArray("P14AHP");
  for (int64_t value : g_3_x) {
    jsonArray.add(value);
  }
  jsonArray = doc.createNestedArray("P15AHP");
  for (int64_t value : h_3_x) {
    jsonArray.add(value);
  }
  jsonArray = doc.createNestedArray("P16AHP");
  jsonArray.add(y_prime);
  jsonArray = doc.createNestedArray("P17AHP");
  jsonArray.add(p_17_AHP);

  jsonArray = doc.createNestedArray("P18AHP");
  for (int64_t value : z) {
    jsonArray.add(value);
  }
  jsonArray = doc.createNestedArray("ComP_AHP_x");
  jsonArray.add(ComP_AHP_x);

  String output;
  serializeJson(doc, output);
  removeFile("/proof.json");
  writeFile("/proof.json", output);
}