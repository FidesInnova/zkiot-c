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

  
  program_param.json
  {
    "Class":  32-bit Integer,
    "A": 64-bit Integer Array,
    "B": 64-bit Integer Array<64-bit Integer Array>
  }

*/

#include "fidesinnova.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "json.hpp"
using ordered_json = nlohmann::ordered_json;
#include <regex>

using namespace std;

int64_t b = 2;
extern "C" void store_register_instances();
// Declare the arrays from assembly as external variables
extern "C" int32_t x0_array[1];
extern "C" int32_t x1_array[1];
extern "C" int32_t x2_array[1];
extern "C" int32_t x3_array[1];
extern "C" int32_t x4_array[1];
extern "C" int32_t x5_array[1];
extern "C" int32_t x6_array[1];
extern "C" int32_t x7_array[1];
extern "C" int32_t x8_array[1];
extern "C" int32_t x9_array[1];
extern "C" int32_t x10_array[1];
extern "C" int32_t x11_array[1];
extern "C" int32_t x12_array[1];
extern "C" int32_t x13_array[1];
extern "C" int32_t x14_array[1];
extern "C" int32_t x15_array[1];
extern "C" int32_t x16_array[1];
extern "C" int32_t x17_array[1];
extern "C" int32_t x18_array[5];
extern "C" int32_t x19_array[1];
extern "C" int32_t x20_array[1];
extern "C" int32_t x21_array[1];
extern "C" int32_t x22_array[1];
extern "C" int32_t x23_array[1];
extern "C" int32_t x24_array[1];
extern "C" int32_t x25_array[1];
extern "C" int32_t x26_array[1];
extern "C" int32_t x27_array[1];
extern "C" int32_t x28_array[1];
extern "C" int32_t x29_array[1];
extern "C" int32_t x30_array[1];
extern "C" int32_t x31_array[1];

void proofGenerator() {
  cout << x18_array[0] << endl;
  // vector<int64_t> z = {1, x0_array[0], x1_array[0], x2_array[0], x3_array[0], x4_array[0], x5_array[0], x6_array[0], x7_array[0], x8_array[0], x9_array[0],
  //                         x10_array[0], x11_array[0], x12_array[0], x13_array[0], x14_array[0], x15_array[0], x16_array[0], x17_array[0], x18_array[0], x19_array[0],
  //                         x20_array[0], x21_array[0], x22_array[0], x23_array[0], x24_array[0], x25_array[0], x26_array[0], x27_array[0], x28_array[0], x29_array[0], 
  //                         x30_array[0], x31_array[0]};
  vector<int64_t> z = {1, 4, 20, 31, 806, 20956};
  // std::ifstream setupFileStream("data/setup3.json");
  // if (!setupFileStream.is_open()) {
  //     std::cerr << "Could not open the setup3.json file!" << std::endl;
  // }
  // nlohmann::json setupJsonData;
  // setupFileStream >> setupJsonData;
  // setupFileStream.close();
  // int64_t Class = setupJsonData["Class"].get<int64_t>();
  // vector<int64_t> ck = setupJsonData["ck"].get<vector<int64_t>>();
  // int64_t vk = setupJsonData["vk"].get<int64_t>();

  // std::ifstream commitmentFileStream("~/zkiot_c_source/gcc/V2.0/data/program_commitment.json");
  // if (!commitmentFileStream.is_open()) {
  //     std::cerr << "Could not open the program_commitment.json file!" << std::endl;
  // }
  // nlohmann::json commitmentJsonData;
  // commitmentFileStream >> commitmentJsonData;
  // commitmentFileStream.close();
  // vector<int64_t> rowA_x = commitmentJsonData["RowA"].get<vector<int64_t>>();
  // vector<int64_t> colA_x = commitmentJsonData["ColA"].get<vector<int64_t>>();
  // vector<int64_t> valA_x = commitmentJsonData["ValA"].get<vector<int64_t>>();
  // vector<int64_t> rowB_x = commitmentJsonData["RowB"].get<vector<int64_t>>();
  // vector<int64_t> colB_x = commitmentJsonData["ColB"].get<vector<int64_t>>();
  // vector<int64_t> valB_x = commitmentJsonData["ValB"].get<vector<int64_t>>();
  // vector<int64_t> rowC_x = commitmentJsonData["RowC"].get<vector<int64_t>>();
  // vector<int64_t> colC_x = commitmentJsonData["ColC"].get<vector<int64_t>>();
  // vector<int64_t> valC_x = commitmentJsonData["ValC"].get<vector<int64_t>>();

  

  // std::ifstream paramFileStream("~/zkiot_c_source/gcc/V2.0/data/program_param.json");
  // if (!paramFileStream.is_open()) {
  //     std::cerr << "Could not open the program_param.json file!" << std::endl;
  // }
  // nlohmann::json paramJsonData;
  // paramFileStream >> paramJsonData;
  // paramFileStream.close();
  // vector<int64_t> nonZeroA = paramJsonData["A"].get<vector<int64_t>>();
  // vector<vector<int64_t>> nonZeroB = paramJsonData["B"].get<vector<vector<int64_t>>>();
  // vector<int64_t> rowA = paramJsonData["rA"].get<vector<int64_t>>();
  // vector<int64_t> colA = paramJsonData["cA"].get<vector<int64_t>>();
  // vector<int64_t> valA = paramJsonData["vA"].get<vector<int64_t>>();
  // vector<int64_t> rowB = paramJsonData["rB"].get<vector<int64_t>>();
  // vector<int64_t> colB = paramJsonData["cB"].get<vector<int64_t>>();
  // vector<int64_t> valB = paramJsonData["vB"].get<vector<int64_t>>();
  // vector<int64_t> rowC = paramJsonData["rC"].get<vector<int64_t>>();
  // vector<int64_t> colC = paramJsonData["cC"].get<vector<int64_t>>();
  // vector<int64_t> valC = paramJsonData["vC"].get<vector<int64_t>>();

  // std::ifstream classFileStream("data/class.json");
  // if (!classFileStream.is_open()) {
  //     std::cerr << "Could not open the file!" << std::endl;
  // }
  // nlohmann::json classJsonData;
  // classFileStream >> classJsonData;
  // classFileStream.close();
  // int64_t n_i, n_g, m, n, p, g;
  // for (const auto& item : classJsonData) {
  //   if (item["Class"] == Class) {
  //     // Number of inputs, gates, m, n, p, and g
  //     n_i = item["n_i"].get<int64_t>();
  //     n_g = item["n_g"].get<int64_t>();
  //     m = item["m"].get<int64_t>();
  //     n = item["n"].get<int64_t>();
  //     p = item["p"].get<int64_t>();
  //     g = item["g"].get<int64_t>();

  //     break; // Stop after finding the first matching "Class"
  //   }
  // }




  int64_t Class = 3;
  vector<int64_t> ck = {11,1309,155771,75218,559337,1106584,774458,1531168,950324,641049,760386,1534921,1396931,81010,1248585,889367,100450,205303,934563,443811,785558,1173747,375250,1018404,350964,1485012,492723,1571123,670006,849627,406353,1363019,1080445,1020559,607409,113868,123724,1296588,1566761,150928,1177222,788775,1556570,616520,1198077,1592199,1499729,565725};
  int64_t vk = 1309;
  vector<int64_t> rowA_x = {469905,454730,902750,1066695,147929,664621,100821,739930};
  vector<int64_t> colA_x = {852571,161650,1543130,388105,25447,115593,41803,1423853};
  vector<int64_t> valA_x = {45582,805822,1515001,845780,376464,555047,619759,350157};
  vector<int64_t> rowB_x = {207638,1648400,628166,1004557,70602,512634,138616,336768};
  vector<int64_t> colB_x = {696691,681295,508950,970277,590324,1035579,1484980,602991};
  vector<int64_t> valB_x = {597622,1323143,17995,122272,1142736,702189,63198,1027496};
  vector<int64_t> rowC_x = {469905,454730,902750,1066695,147929,664621,100821,739930};
  vector<int64_t> colC_x = {469905,454730,902750,1066695,147929,664621,100821,739930};
  vector<int64_t> valC_x = {1083027,544916,467855,1126753,525698,1158417,1188860,435354};

  vector<int64_t> nonZeroA = {19,0,19,19};
  vector<vector<int64_t>> nonZeroB = {{33,20,1},{34,0,11},{34,19,1},{35,21,1},{36,21,1}};
  vector<int64_t> rowA = {1190739,1632256,159545,899305,373750,668759,747302,1444226};
  vector<int64_t> colA = {1195510,1,1195510,1195510,373750,668759,747302,1444226};
  vector<int64_t> valA = {78649,1088607,1609535,944507,0,0,0,0};
  vector<int64_t> rowB = {1190739,1632256,1632256,159545,899305,668759,747302,1444226};
  vector<int64_t> colB = {1536124,1,1195510,1669124,1669124,668759,747302,1444226};
  vector<int64_t> valB = {1640009,226430,1640009,949756,324772,0,0,0};
  vector<int64_t> rowC = {1190739,1632256,159545,899305,373750,668759,747302,1444226};
  vector<int64_t> colC = {1190739,1632256,159545,899305,373750,668759,747302,1444226};
  vector<int64_t> valC = {1495917,1550025,1582341,679291,0,0,0,0};

  int64_t n_i, n_g, m, n, p, g;
  n_i = 32;
  n_g = 4;
  m = 8;
  n = 37;
  p = 1678321;
  g = 11;



  int64_t t = n_i + 1;


  // Initialize matrices A, B, C
  vector<vector<int64_t>> A(n, vector<int64_t>(n, 0ll));
  vector<vector<int64_t>> B(n, vector<int64_t>(n, 0ll));
  vector<vector<int64_t>> C(n, vector<int64_t>(n, 0ll));

  int64_t rowMatA = n_i;
  for (int64_t i = 0; i < nonZeroA.size(); i++) {
    int64_t col = nonZeroA[i];
    // Set the value in the matrix
    A[i][col] = 1;
  }
  for (int64_t i = 0; i < nonZeroB.size(); i++) {
    int64_t row = nonZeroB[i][0];
    int64_t col = nonZeroB[i][1];
    int64_t val = nonZeroB[i][2];
    
    // Set the value in the matrix
    B[row][col] = val;
  }
  for(int64_t i = (n - n_g); i < n; i++) {
    // Set the value in the matrix
    C[i][i] = 1;
  }
  Polynomial::printMatrix(A, "A");
  Polynomial::printMatrix(B, "B");
  Polynomial::printMatrix(C, "C");

  // vector<int64_t> z;

  vector<int64_t> H;
  int64_t w, g_n;

  H.push_back(1);
  g_n = ((p - 1) / n) % p;
  w = Polynomial::power(g, g_n, p);
  for (int64_t i = 1; i < n; i++) {
    H.push_back(Polynomial::power(w, i, p));
  }
  cout << "H[n]: ";
  for (int64_t i = 0; i < n; i++) {
    cout << H[i] << " ";
  }
  cout << endl;
  
  int64_t y, g_m;

  vector<int64_t> K;
  K.push_back(1);
  g_m = ((p - 1) * Polynomial::pInverse(m, p)) % p;
  y = Polynomial::power(g, g_m, p);
  for (int64_t i = 1; i < m; i++) {
    K.push_back(Polynomial::power(y, i, p));
  }
  cout << "K[m]: ";
  for (int64_t i = 0; i < m; i++) {
    cout << K[i] << " ";
  }
  cout << endl;


  vector<vector<int64_t>> Az(n, vector<int64_t>(1, 0));
  vector<vector<int64_t>> Bz(n, vector<int64_t>(1, 0));
  vector<vector<int64_t>> Cz(n, vector<int64_t>(1, 0));
  cout << "n_i: " << n_i << endl;
  
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
  cout << "Matrice Az under modulo " << p << " is: ";
  for (int64_t i = 0; i < n; i++) {
    cout << Az[i][0] << " ";
  }
  cout << endl;

  cout << "Matrice Bz under modulo " << p << " is: ";
  for (int64_t i = 0; i < n; i++) {
    cout << Bz[i][0] << " ";
  }
  cout << endl;

  cout << "Matrice Cz under modulo " << p << " is: ";
  for (int64_t i = 0; i < n; i++) {
    cout << Cz[i][0] << " ";
  }
  cout << endl;



  vector<vector<int64_t>> zA(2);
  cout << "zA(x):" << endl;
  for (int64_t i = 0; i < n + b; i++) {
  // for (int64_t i = 0; i < n; i++) {
    if (i < n) {
      zA[0].push_back(H[i]);
      zA[1].push_back(Az[i][0]);
    } else {
      zA[0].push_back(Polynomial::generateRandomNumber(H, p - n));
      zA[1].push_back(Polynomial::generateRandomNumber(H, p - n));
    }
    cout << "zA(" << zA[0][i] << ")= " << zA[1][i];
  }
  // zA[0].push_back(150);
  // zA[1].push_back(5);
  // zA[0].push_back(80);
  // zA[1].push_back(47);

  vector<int64_t> z_hatA = Polynomial::setupLagrangePolynomial(zA[0], zA[1], p, "z_hatA(x)");


  vector<vector<int64_t>> zB(2);
  cout << "zB(x):" << endl;
  for (int64_t i = 0; i < n + b; i++) {
  // for (int64_t i = 0; i < n; i++) {
    if (i < n) {
      zB[0].push_back(H[i]);
      zB[1].push_back(Bz[i][0]);
    } else {
      zB[0].push_back(zA[0][i]);
      zB[1].push_back(Polynomial::generateRandomNumber(H, p - n));
    }
    cout << "zB(" << zB[0][i] << ")= " << zB[1][i] << endl;
  }
  // zB[0].push_back(150);
  // zB[1].push_back(15);
  // zB[0].push_back(80);
  // zB[1].push_back(170);
  vector<int64_t> z_hatB = Polynomial::setupLagrangePolynomial(zB[0], zB[1], p, "z_hatB(x)");

  vector<vector<int64_t>> zC(2);
  cout << "zC(x):";
  for (int64_t i = 0; i < n + b; i++) {
  // for (int64_t i = 0; i < n; i++) {
    if (i < n) {
      zC[0].push_back(H[i]);
      zC[1].push_back(Cz[i][0]);
    } else {
      zC[0].push_back(zA[0][i]);
      zC[1].push_back(Polynomial::generateRandomNumber(H, p - n));
    }
    cout << "zC(" << zC[0][i] << ")= " << zC[1][i] << endl;
  }
  // zC[0].push_back(150);
  // zC[1].push_back(1);
  // zC[0].push_back(80);
  // zC[1].push_back(100);
  vector<int64_t> z_hatC = Polynomial::setupLagrangePolynomial(zC[0], zC[1], p, "z_hatC(x)");


  vector<int64_t> zero_to_t_for_H;
  vector<int64_t> t_to_n_for_H;
  vector<int64_t> zero_to_t_for_z;
  vector<int64_t> t_to_n_for_z;
  for (int64_t i = 0; i < t; i++) {
    zero_to_t_for_H.push_back(H[i]);
    zero_to_t_for_z.push_back(z[i]);
  }
  vector<int64_t> polyX_HAT_H = Polynomial::setupLagrangePolynomial(zero_to_t_for_H, zero_to_t_for_z, p, "x_hat(h)");

  for (int64_t i = 0; i < n - t; i++) {
    t_to_n_for_H.push_back(H[i + t]);
    t_to_n_for_z.push_back(z[i + t]);
  }
  cout << "w_bar(h):" << endl;
  vector<int64_t> w_bar(n - t + b);
  vector<int64_t> w_bar_numerator(n - t, 1);
  vector<int64_t> w_bar_denominator(n - t, 1);
  for (int64_t i = 0; i < n - t; i++) {
    w_bar_numerator[i] = (t_to_n_for_z[i] - (Polynomial::evaluatePolynomial(polyX_HAT_H, t_to_n_for_H[i], p)));

    for (int64_t j = 0; j < zero_to_t_for_H.size(); j++) {
      w_bar_denominator[i] *= (t_to_n_for_H[i] - zero_to_t_for_H[j]);
      // Apply pulus to keep the number within the bounds
      w_bar_denominator[i] %= p;

      // // If result becomes negative, convert it to a positive equivalent under modulo
      if (w_bar_denominator[i] < 0) {
        w_bar_denominator[i] += p;
      }
    }
    w_bar_denominator[i] = Polynomial::pInverse(w_bar_denominator[i], p);
    w_bar[i] = (w_bar_numerator[i] * w_bar_denominator[i]) % p;

    if (w_bar[i] < 0) {
      w_bar[i] += p;
    }
    cout << "w_bar(" << t_to_n_for_H[i] << ")= " << w_bar[i] << endl;
  }

  cout << "w_hat(x):" << endl;
  vector<vector<int64_t>> w_hat(2);
  for (int64_t i = 0; i < n - t; i++) {
    w_hat[0].push_back(t_to_n_for_H[i]);
    w_hat[1].push_back(w_bar[i]);
    cout << "w_hat(" << w_hat[0][i] << ")= " << w_hat[1][i] << endl;
  }
  // w_hat[0].push_back(150);
  // w_hat[1].push_back(42);
  // w_hat[0].push_back(80);
  // w_hat[1].push_back(180);
  for (int64_t i = n; i < n + b; i++) {
    w_hat[0].push_back(zA[0][i]);
    w_hat[1].push_back(Polynomial::generateRandomNumber(H, p));
    cout << "w_hat(" << w_hat[0][i - b] << ")= " << w_hat[1][i - b] << endl;
  }
  vector<int64_t> w_hat_x = Polynomial::setupLagrangePolynomial(w_hat[0], w_hat[1], p, "w_hat(x)");

  vector<int64_t> productAB = Polynomial::multiplyPolynomials(z_hatA, z_hatB, p);
  vector<int64_t> zAzB_zC = Polynomial::subtractPolynomials(productAB, z_hatC, p);
  Polynomial::printPolynomial(zAzB_zC, "zA(x)zB(x)-zC(x)");


  vector<int64_t> vH_x(n + 1, 0);
  vH_x[0] = (-1) % p;
  if (vH_x[0] < 0) {
    vH_x[0] += p; //Handling negative values properly
  }
  vH_x[n] = 1;
  Polynomial::printPolynomial(vH_x, "KvH(x)");


  vector<int64_t> vK_x = Polynomial::createLinearPolynomial(K[0]);
  // Multiply (x - K) for all other K
  for (size_t i = 1; i < K.size(); i++) {
    vector<int64_t> nextPoly = Polynomial::createLinearPolynomial(K[i]);
    vK_x = Polynomial::multiplyPolynomials(vK_x, nextPoly, p);
  }
  Polynomial::printPolynomial(vK_x, "vK(x)");

  // Dividing the product of zAzB_zC by vH_x
  vector<int64_t> h_0_x = Polynomial::dividePolynomials(zAzB_zC, vH_x, p)[0];
  Polynomial::printPolynomial(h_0_x, "h0(x)");

  // vector<int64_t> s_x = generateRandomPolynomial(n, (2*n)+b-1, p);
  vector<int64_t> s_x = { 115, 3, 0, 0, 20, 1, 0, 17, 101, 0, 5 };
  Polynomial::printPolynomial(s_x, "s(x)");

  int64_t sigma1 = Polynomial::sumOfEvaluations(s_x, H, p);
  cout << "sigma1 = " << sigma1 << endl;

  int64_t alpha = 10;
  int64_t beta1 = 22;  // int64_t beta1 = hashAndExtractLower4Bytes(s_x[8], p);
  int64_t beta2 = 80;  // int64_t beta2 = hashAndExtractLower4Bytes(s_x[9], p);
  int64_t etaA = 2;
  int64_t etaB = 30;
  int64_t etaC = 100;

  cout << "alpha = " << alpha << endl;
  cout << "beta1 = " << beta1 << endl;
  cout << "beta2 = " << beta2 << endl;
  cout << "etaA = " << etaA << endl;
  cout << "etaB = " << etaB << endl;
  cout << "etaC = " << etaC << endl;


  vector<int64_t> etaA_z_hatA_x = Polynomial::multiplyPolynomialByNumber(z_hatA, etaA, p);
  vector<int64_t> etaB_z_hatB_x = Polynomial::multiplyPolynomialByNumber(z_hatB, etaB, p);
  vector<int64_t> etaC_z_hatC_x = Polynomial::multiplyPolynomialByNumber(z_hatC, etaC, p);
  Polynomial::printPolynomial(etaA_z_hatA_x, "etaA_z_hatA(x)");
  Polynomial::printPolynomial(etaB_z_hatB_x, "etaB_z_hatB(x)");
  Polynomial::printPolynomial(etaC_z_hatC_x, "etaC_z_hatC(x)");

  vector<int64_t> Sum_M_eta_M_z_hat_M_x = Polynomial::addPolynomials(Polynomial::addPolynomials(etaA_z_hatA_x, etaB_z_hatB_x, p), etaC_z_hatC_x, p);
  Polynomial::printPolynomial(Sum_M_eta_M_z_hat_M_x, "Sum_M_z_hatM(x)");

  vector<int64_t> r_alpha_x = Polynomial::calculatePolynomial_r_alpha_x(alpha, n, p);
  Polynomial::printPolynomial(r_alpha_x, "r(alpha, x)");

  vector<int64_t> r_Sum_x = Polynomial::multiplyPolynomials(r_alpha_x, Sum_M_eta_M_z_hat_M_x, p);
  Polynomial::printPolynomial(r_Sum_x, "r(alpha, x)Sum_M_z_hatM(x)");

  vector<int64_t> v_H = Polynomial::expandPolynomials(zero_to_t_for_H, p);
  Polynomial::printPolynomial(v_H, "v_H");
  vector<int64_t> z_hat_x = Polynomial::addPolynomials(Polynomial::multiplyPolynomials(w_hat_x, v_H, p), polyX_HAT_H, p);
  Polynomial::printPolynomial(z_hat_x, "z_hat(x)");

  vector<int64_t> A_hat(2);
  for (int64_t i = 0; i < nonZeroA.size(); i++) {
    vector<int64_t> buff_n = Polynomial::calculatePolynomial_r_alpha_x(rowA[i], H.size(), p);
    int64_t eval = Polynomial::evaluatePolynomial(buff_n, rowA[i], p);
    vector<int64_t> buff = Polynomial::calculatePolynomial_r_alpha_x(colA[i], H.size(), p);
    eval = (eval * valA[i]) % p;
    if (eval < 0) {
      eval += p;
    }
    buff = Polynomial::multiplyPolynomialByNumber(buff, eval, p);
    buff = Polynomial::multiplyPolynomialByNumber(buff, Polynomial::calculatePolynomial_r_alpha_k(alpha, rowA[i], H.size(), p), p);
    if (i > 0) {
      A_hat = Polynomial::addPolynomials(A_hat, buff, p);
    } else {
      A_hat = buff;
    }
  }
  Polynomial::printPolynomial(A_hat, "A_hat(x)");

  vector<int64_t> B_hat(2);
  for (int64_t i = 0; i < nonZeroB[0].size(); i++) {
    vector<int64_t> buff_n = Polynomial::calculatePolynomial_r_alpha_x(rowB[i], H.size(), p);
    int64_t eval = Polynomial::evaluatePolynomial(buff_n, rowB[i], p);
    vector<int64_t> buff = Polynomial::calculatePolynomial_r_alpha_x(colB[i], H.size(), p);
    eval = (eval * valB[i]) % p;
    if (eval < 0) {
      eval += p;
    }
    buff = Polynomial::multiplyPolynomialByNumber(buff, eval, p);
    buff = Polynomial::multiplyPolynomialByNumber(buff, Polynomial::calculatePolynomial_r_alpha_k(alpha, rowB[i], H.size(), p), p);
    if (i > 0) {
      B_hat = Polynomial::addPolynomials(B_hat, buff, p);
    } else {
      B_hat = buff;
    }
  }
  Polynomial::printPolynomial(B_hat, "B_hat(x)");
  
  vector<int64_t> C_hat(2);
  for (int64_t i = 0; i < n_g; i++) {
    vector<int64_t> buff_n = Polynomial::calculatePolynomial_r_alpha_x(rowC[i], H.size(), p);
    int64_t eval = Polynomial::evaluatePolynomial(buff_n, rowC[i], p);
    vector<int64_t> buff = Polynomial::calculatePolynomial_r_alpha_x(colC[i], H.size(), p);
    eval = (eval * valC[i]) % p;
    if (eval < 0) {
      eval += p;
    }
    buff = Polynomial::multiplyPolynomialByNumber(buff, eval, p);
    buff = Polynomial::multiplyPolynomialByNumber(buff, Polynomial::calculatePolynomial_r_alpha_k(alpha, rowC[i], H.size(), p), p);
    if (i > 0) {
      C_hat = Polynomial::addPolynomials(C_hat, buff, p);
    } else {
      C_hat = buff;
    }
  }
  Polynomial::printPolynomial(C_hat, "C_hat(x)");

/*
// Check for overflow risk for etaA, etaB and etaC and ensure p is applied correctly
// vector<int64_t> etaA_z_hatA_x = Polynomial::multiplyPolynomialByNumber(A_hatA, etaA % p, p);
// vector<int64_t> etaB_z_hatB_x = Polynomial::multiplyPolynomialByNumber(B_hatB, etaB % p, p);
// vector<int64_t> etaC_z_hatC_x = Polynomial::multiplyPolynomialByNumber(C_hatC, etaC % p, p);*/

  vector<int64_t> eta_A_hat = Polynomial::multiplyPolynomialByNumber(A_hat, etaA, p);
  vector<int64_t> eta_B_hat = Polynomial::multiplyPolynomialByNumber(B_hat, etaB, p);
  vector<int64_t> eta_C_hat = Polynomial::multiplyPolynomialByNumber(C_hat, etaC, p);
  Polynomial::printPolynomial(eta_A_hat, "eta_A_hat: ");
  Polynomial::printPolynomial(eta_B_hat, "eta_B_hat: ");
  Polynomial::printPolynomial(eta_C_hat, "eta_C_hat: ");

  // Calculate the sum of the three polynomials and print the result
  vector<int64_t> Sum_M_eta_M_r_M_alpha_x = Polynomial::addPolynomials(Polynomial::addPolynomials(eta_A_hat, eta_B_hat, p), eta_C_hat, p);
  Polynomial::printPolynomial(Sum_M_eta_M_r_M_alpha_x, "Sum_M_eta_M_r_M(alpha ,x)");

  // Multiply the sum by another polynomial z_hat_x and print the result
  vector<int64_t> Sum_M_eta_M_r_M_alpha_x_z_hat_x = Polynomial::multiplyPolynomials(Sum_M_eta_M_r_M_alpha_x, z_hat_x, p);
  Polynomial::printPolynomial(Sum_M_eta_M_r_M_alpha_x_z_hat_x, "Sum_M_eta_M_r_M(alpha ,x)z-hat(x)");

  // Calculate the sum for the check protocol, subtracting the pified sum from s_x
  vector<int64_t> Sum_check_protocol = Polynomial::addPolynomials(s_x, (Polynomial::subtractPolynomials(r_Sum_x, Sum_M_eta_M_r_M_alpha_x_z_hat_x, p)), p);
  Polynomial::printPolynomial(Sum_check_protocol, "Sum_check_protocol");

  // Divide the sum check protocol by vH_x to get two results: h1(x) and g1(x)
  vector<int64_t> h_1_x = Polynomial::dividePolynomials(Sum_check_protocol, vH_x, p)[0];
  Polynomial::printPolynomial(h_1_x, "h1(x)");

  // Get the second part of the division result, g1(x), and erase the first element
  vector<int64_t> g_1_x = Polynomial::dividePolynomials(Sum_check_protocol, vH_x, p)[1];
  g_1_x.erase(g_1_x.begin());
  Polynomial::printPolynomial(g_1_x, "g1(x)");

  // Calculate sigma2 using the evaluations of the polynomials A_hat, B_hat, and C_hat and print the result of sigma2
  int64_t sigma2 = ((etaA * Polynomial::evaluatePolynomial(A_hat, beta1, p)) % p + (etaB * Polynomial::evaluatePolynomial(B_hat, beta1, p)) % p + (etaC * Polynomial::evaluatePolynomial(C_hat, beta1, p)) % p) % p;
  cout << "sigma2 = " << sigma2 << endl;

  // Initialize vectors for the pified polynomial results with zeros
  vector<int64_t> A_hat_M_hat(H.size(), 0);
  vector<int64_t> B_hat_M_hat(H.size(), 0);
  vector<int64_t> C_hat_M_hat(H.size(), 0);

  // Loop through non-zero rows for matrix A and calculate the pified polynomial A_hat_M_hat
  for (int64_t i = 0; i < nonZeroA.size(); i++) {
    vector<int64_t> buff_nA = Polynomial::calculatePolynomial_r_alpha_x(colA[i], H.size(), p);
    int64_t evalA = Polynomial::evaluatePolynomial(buff_nA, beta1, p);
    vector<int64_t> buffA = Polynomial::calculatePolynomial_r_alpha_x(rowA[i], H.size(), p);
    evalA = (evalA * valA[i]) % p;
    if (evalA < 0) {
      evalA += p;
    }
    buffA = Polynomial::multiplyPolynomialByNumber(buffA, evalA, p);

    if (i > 0) {
      A_hat_M_hat = Polynomial::addPolynomials(A_hat_M_hat, buffA, p);
    } else {
      A_hat_M_hat = buffA;
    }
  }
  for (int64_t i = 0; i < nonZeroB[0].size(); i++) {
    vector<int64_t> buff_nB = Polynomial::calculatePolynomial_r_alpha_x(colB[i], H.size(), p);
    int64_t evalB = Polynomial::evaluatePolynomial(buff_nB, beta1, p);
    vector<int64_t> buffB = Polynomial::calculatePolynomial_r_alpha_x(rowB[i], H.size(), p);
    evalB = (evalB * valB[i]) % p;
    if (evalB < 0) {
      evalB += p;
    }
    buffB = Polynomial::multiplyPolynomialByNumber(buffB, evalB, p);

    if (i > 0) {
      B_hat_M_hat = Polynomial::addPolynomials(B_hat_M_hat, buffB, p);
    } else {
      B_hat_M_hat = buffB;
    }
  }
  for (int64_t i = 0; i < n_g; i++) {
    vector<int64_t> buff_nC = Polynomial::calculatePolynomial_r_alpha_x(colC[i], H.size(), p);
    int64_t evalC = Polynomial::evaluatePolynomial(buff_nC, beta1, p);
    vector<int64_t> buffC = Polynomial::calculatePolynomial_r_alpha_x(rowC[i], H.size(), p);
    evalC = (evalC * valC[i]) % p;
    if (evalC < 0) {
      evalC += p;
    }
    buffC = Polynomial::multiplyPolynomialByNumber(buffC, evalC, p);

    if (i > 0) {
      C_hat_M_hat = Polynomial::addPolynomials(C_hat_M_hat, buffC, p);
    } else {
      C_hat_M_hat = buffC;
    }
  }
  // Print the final pified polynomials for A, B, and C
  Polynomial::printPolynomial(A_hat_M_hat, "A_hat_M_hat");
  Polynomial::printPolynomial(B_hat_M_hat, "B_hat_M_hat");
  Polynomial::printPolynomial(C_hat_M_hat, "C_hat_M_hat");

  // Multiply the pified polynomials by their respective eta values and print
  vector<int64_t> eta_A_hat_M_hat = Polynomial::multiplyPolynomialByNumber(A_hat_M_hat, etaA, p);
  vector<int64_t> eta_B_hat_M_hat = Polynomial::multiplyPolynomialByNumber(B_hat_M_hat, etaB, p);
  vector<int64_t> eta_C_hat_M_hat = Polynomial::multiplyPolynomialByNumber(C_hat_M_hat, etaC, p);
  Polynomial::printPolynomial(eta_A_hat_M_hat, "eta_A_hat_M_hat: ");
  Polynomial::printPolynomial(eta_B_hat_M_hat, "eta_B_hat_M_hat: ");
  Polynomial::printPolynomial(eta_C_hat_M_hat, "eta_C_hat_M_hat: ");

  // Calculate the final result for r_Sum_M_eta_M_M_hat_x_beta1
  vector<int64_t> r_Sum_M_eta_M_M_hat_x_beta1 = Polynomial::multiplyPolynomials(Polynomial::addPolynomials(Polynomial::addPolynomials(eta_A_hat_M_hat, eta_B_hat_M_hat, p), eta_C_hat_M_hat, p), r_alpha_x, p);
  Polynomial::printPolynomial(r_Sum_M_eta_M_M_hat_x_beta1, "r_Sum_M_eta_M_M_hat_x_beta1");

  // Divide the final result by vH_x to get h2(x) and g2(x)
  vector<int64_t> h_2_x = Polynomial::dividePolynomials(r_Sum_M_eta_M_M_hat_x_beta1, vH_x, p)[0];
  Polynomial::printPolynomial(h_2_x, "h2(x)");

  vector<int64_t> g_2_x = Polynomial::dividePolynomials(r_Sum_M_eta_M_M_hat_x_beta1, vH_x, p)[1];
  g_2_x.erase(g_2_x.begin());//remove the first item
  Polynomial::printPolynomial(g_2_x, "g2(x)");

  // Generate Lagrange polynomials for row, column, and value vectors for matrices A and B
  // vector<int64_t> rowA_x = Polynomial::setupLagrangePolynomial(rowA[0], rowA[1], p, "rowA(x)");
  // vector<int64_t> colA_x = Polynomial::setupLagrangePolynomial(colA[0], colA[1], p, "colA(x)");
  // vector<int64_t> valA_x = Polynomial::setupLagrangePolynomial(valA[0], valA[1], p, "valA(x)");

  // vector<int64_t> rowB_x = Polynomial::setupLagrangePolynomial(rowB[0], rowB[1], p, "rowB(x)");
  // vector<int64_t> colB_x = Polynomial::setupLagrangePolynomial(colB[0], colB[1], p, "colB(x)");
  // vector<int64_t> valB_x = Polynomial::setupLagrangePolynomial(valB[0], valB[1], p, "valB(x)");

  // vector<int64_t> rowC_x = Polynomial::setupLagrangePolynomial(rowC[0], rowC[1], p, "rowC(x)");
  // vector<int64_t> colC_x = Polynomial::setupLagrangePolynomial(colC[0], colC[1], p, "colC(x)");
  // vector<int64_t> valC_x = Polynomial::setupLagrangePolynomial(valC[0], valC[1], p, "valC(x)");

  // Evaluate polynomial vH at beta1 and beta2
  int64_t vH_beta1 = Polynomial::evaluatePolynomial(vH_x, beta1, p);
  cout << "vH(beta1) = " << vH_beta1 << endl;

  int64_t vH_beta2 = Polynomial::evaluatePolynomial(vH_x, beta2, p);
  cout << "vH(beta2) = " << vH_beta2 << endl;

  // Initialize vectors for function points and sigma value
  vector<int64_t> points_f_3(K.size(), 0);
  int64_t sigma3 = 0;

  // Loop over K to compute delta and signature values for A, B, and C
  for (int64_t i = 0; i < K.size(); i++) {
    int64_t deA = (beta2 - Polynomial::evaluatePolynomial(rowA_x, K[i], p)) * (beta1 - Polynomial::evaluatePolynomial(colA_x, K[i], p)) % p;
    if (deA < 0) deA += p;

    int64_t deB = (beta2 - Polynomial::evaluatePolynomial(rowB_x, K[i], p)) * (beta1 - Polynomial::evaluatePolynomial(colB_x, K[i], p)) % p;
    if (deB < 0) deB += p;

    int64_t deC = (beta2 - Polynomial::evaluatePolynomial(rowC_x, K[i], p)) * (beta1 - Polynomial::evaluatePolynomial(colC_x, K[i], p)) % p;
    if (deC < 0) deC += p;

    int64_t sig3_A = (etaA * (vH_beta2 * vH_beta1 % p) % p * Polynomial::evaluatePolynomial(valA_x, K[i], p) % p * Polynomial::pInverse(deA, p)) % p;
    if (sig3_A < 0) sig3_A += p;

    int64_t sig3_B = (etaB * (vH_beta2 * vH_beta1 % p) % p * Polynomial::evaluatePolynomial(valB_x, K[i], p) % p * Polynomial::pInverse(deB, p)) % p;
    if (sig3_B < 0) sig3_B += p;

    int64_t sig3_C = (etaC * (vH_beta2 * vH_beta1 % p) % p * Polynomial::evaluatePolynomial(valC_x, K[i], p) % p * Polynomial::pInverse(deC, p)) % p;
    if (sig3_C < 0) sig3_C += p;

    points_f_3[i] = (sig3_A + sig3_B + sig3_C) % p;
    sigma3 += points_f_3[i];
    sigma3 %= p;
  }
  cout << "sigma3 = " << sigma3 << endl;

  // Create polynomials for beta1 and beta2
  vector<int64_t> poly_beta1 = { beta1 };
  vector<int64_t> poly_beta2 = { beta2 };

  // Compute polynomial products for signatures
  vector<int64_t> poly_pi_a = Polynomial::multiplyPolynomials(Polynomial::subtractPolynomials(rowA_x, poly_beta2, p), Polynomial::subtractPolynomials(colA_x, poly_beta1, p), p);
  vector<int64_t> poly_pi_b = Polynomial::multiplyPolynomials(Polynomial::subtractPolynomials(rowB_x, poly_beta2, p), Polynomial::subtractPolynomials(colB_x, poly_beta1, p), p);
  vector<int64_t> poly_pi_c = Polynomial::multiplyPolynomials(Polynomial::subtractPolynomials(rowC_x, poly_beta2, p), Polynomial::subtractPolynomials(colC_x, poly_beta1, p), p);
  Polynomial::printPolynomial(poly_pi_c, "poly_pi_a");
  Polynomial::printPolynomial(poly_pi_b, "poly_pi_b");
  Polynomial::printPolynomial(poly_pi_c, "poly_pi_c");

  // Compute polynomials for signature multipliers
  vector<int64_t> poly_etaA_vH_B2_vH_B1 = { (etaA * vH_beta2 * vH_beta1) % p };
  vector<int64_t> poly_etaB_vH_B2_vH_B1 = { (etaB * vH_beta2 * vH_beta1) % p };
  vector<int64_t> poly_etaC_vH_B2_vH_B1 = { (etaC * vH_beta2 * vH_beta1) % p };

  // Calculate signatures
  vector<int64_t> poly_sig_a = Polynomial::multiplyPolynomials(poly_etaA_vH_B2_vH_B1, valA_x, p);
  vector<int64_t> poly_sig_b = Polynomial::multiplyPolynomials(poly_etaB_vH_B2_vH_B1, valB_x, p);
  vector<int64_t> poly_sig_c = Polynomial::multiplyPolynomials(poly_etaC_vH_B2_vH_B1, valC_x, p);

  vector<int64_t> a_x = Polynomial::addPolynomials(Polynomial::addPolynomials(Polynomial::multiplyPolynomials(poly_sig_a, Polynomial::multiplyPolynomials(poly_pi_b, poly_pi_c, p), p), Polynomial::multiplyPolynomials(poly_sig_b, Polynomial::multiplyPolynomials(poly_pi_a, poly_pi_c, p), p), p), Polynomial::multiplyPolynomials(poly_sig_c, Polynomial::multiplyPolynomials(poly_pi_a, poly_pi_b, p), p), p);
  Polynomial::printPolynomial(a_x, "a(x)");

  vector<int64_t> b_x = Polynomial::multiplyPolynomials(Polynomial::multiplyPolynomials(poly_pi_a, poly_pi_b, p), poly_pi_c, p);
  Polynomial::printPolynomial(b_x, "b(x)");

  // Set up polynomial for f_3 using K
  vector<int64_t> poly_f_3x = Polynomial::setupLagrangePolynomial(K, points_f_3, p, "poly_f_3(x)");

  vector<int64_t> g_3_x = poly_f_3x;
  g_3_x.erase(g_3_x.begin());
  Polynomial::printPolynomial(g_3_x, "g3(x)");

  // Calculate sigma_3_set_k based on sigma3 and K.size()
  vector<int64_t> sigma_3_set_k;
  sigma_3_set_k.push_back((sigma3 * Polynomial::pInverse(K.size(), p)) % p);
  cout << "sigma_3_set_k = " << sigma_3_set_k[0] << endl;

  // Update polynomial f_3 by subtracting sigma_3_set_k
  vector<int64_t> poly_f_3x_new = Polynomial::subtractPolynomials(poly_f_3x, sigma_3_set_k, p);
  Polynomial::printPolynomial(poly_f_3x_new, "f3(x)new");

  // Calculate polynomial h_3(x) using previous results
  vector<int64_t> h_3_x = Polynomial::dividePolynomials(Polynomial::subtractPolynomials(a_x, Polynomial::multiplyPolynomials(b_x, Polynomial::addPolynomials(poly_f_3x_new, sigma_3_set_k, p), p), p), vK_x, p)[0];
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
      Polynomial::addPolynomials(Polynomial::addPolynomials(Polynomial::multiplyPolynomialByNumber(w_hat_x, eta_w_hat, p), Polynomial::multiplyPolynomialByNumber(z_hatA, eta_z_hatA, p), p),
                                 Polynomial::addPolynomials(Polynomial::multiplyPolynomialByNumber(z_hatB, eta_z_hatB, p), Polynomial::multiplyPolynomialByNumber(z_hatC, eta_z_hatC, p), p), p),
      Polynomial::addPolynomials(
        Polynomial::addPolynomials(Polynomial::addPolynomials(Polynomial::multiplyPolynomialByNumber(h_0_x, eta_h_0_x, p), Polynomial::multiplyPolynomialByNumber(s_x, eta_s_x, p), p),
                                   Polynomial::addPolynomials(Polynomial::multiplyPolynomialByNumber(g_1_x, eta_g_1_x, p), Polynomial::multiplyPolynomialByNumber(h_1_x, eta_h_1_x, p), p), p),
        Polynomial::addPolynomials(Polynomial::addPolynomials(Polynomial::multiplyPolynomialByNumber(g_2_x, eta_g_2_x, p), Polynomial::multiplyPolynomialByNumber(h_2_x, eta_h_2_x, p), p),
                                   Polynomial::addPolynomials(Polynomial::multiplyPolynomialByNumber(g_3_x, eta_g_3_x, p), Polynomial::multiplyPolynomialByNumber(h_3_x, eta_h_3_x, p), p), p),
        p),
      p);
  Polynomial::printPolynomial(p_x, "p(x)");

  // Choose a random value for z
  int64_t z_random = 2;
  int64_t y_prime = Polynomial::evaluatePolynomial(p_x, z_random, p);
  cout << "y_prime = " << y_prime << endl;

  // p(x) - y'  =>  p(x)  !!!!!!!
  vector<int64_t> q_xBuf;
  q_xBuf.push_back(p - z_random);
  q_xBuf.push_back(1);
  Polynomial::printPolynomial(q_xBuf, "div = ");

  vector<int64_t> q_x = Polynomial::dividePolynomials(p_x, q_xBuf, p)[0];
  Polynomial::printPolynomial(q_x, "q(x)");

  // Generate a KZG commitment for q(x) using the provided verification key (ck)
  int64_t p_17_AHP = Polynomial::KZG_Commitment(ck, q_x, p);
  cout << "p_17_AHP = " << p_17_AHP << endl;

  // Generate KZG commitments for various polynomials and print
  int64_t Com1_AHP_x = Polynomial::KZG_Commitment(ck, w_hat_x, p);
  int64_t Com2_AHP_x = Polynomial::KZG_Commitment(ck, z_hatA, p);
  int64_t Com3_AHP_x = Polynomial::KZG_Commitment(ck, z_hatB, p);
  int64_t Com4_AHP_x = Polynomial::KZG_Commitment(ck, z_hatC, p);
  int64_t Com5_AHP_x = Polynomial::KZG_Commitment(ck, h_0_x, p);
  int64_t Com6_AHP_x = Polynomial::KZG_Commitment(ck, s_x, p);
  int64_t Com7_AHP_x = Polynomial::KZG_Commitment(ck, g_1_x, p);
  int64_t Com8_AHP_x = Polynomial::KZG_Commitment(ck, h_1_x, p);
  int64_t Com9_AHP_x = Polynomial::KZG_Commitment(ck, g_2_x, p);
  int64_t Com10_AHP_x = Polynomial::KZG_Commitment(ck, h_2_x, p);
  int64_t Com11_AHP_x = Polynomial::KZG_Commitment(ck, g_3_x, p);
  int64_t Com12_AHP_x = Polynomial::KZG_Commitment(ck, h_3_x, p);
  cout << "Com1_AHP_x = " << Com1_AHP_x << endl;
  cout << "Com2_AHP_x = " << Com2_AHP_x << endl;
  cout << "Com3_AHP_x = " << Com3_AHP_x << endl;
  cout << "Com4_AHP_x = " << Com4_AHP_x << endl;
  cout << "Com5_AHP_x = " << Com5_AHP_x << endl;
  cout << "Com6_AHP_x = " << Com6_AHP_x << endl;
  cout << "Com7_AHP_x = " << Com7_AHP_x << endl;
  cout << "Com8_AHP_x = " << Com8_AHP_x << endl;
  cout << "Com9_AHP_x = " << Com9_AHP_x << endl;
  cout << "Com10_AHP_x = " << Com10_AHP_x << endl;
  cout << "Com11_AHP_x = " << Com11_AHP_x << endl;
  cout << "Com12_AHP_x = " << Com12_AHP_x << endl;

  // Generate a KZG commitment for the combined polynomial p(x)
  int64_t ComP_AHP_x = Polynomial::KZG_Commitment(ck, p_x, p);
  cout << "ComP_AHP = " << ComP_AHP_x << endl;



  // TODO: Add req params to the proof.json
  /*
    "CommitmentID":  64-bit,
    "Class": 32-bit Integer,
    "IoT_Manufacturer_Name": String,
    "IoT_Device_Name": String,
    "Device_Hardware_Version": float,
    "Firmware_Version": float,
    "code_block": 64-bit Array,
  */
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
  proof["P16AHP"] = y_prime;
  proof["P17AHP"] = p_17_AHP;

  cout << "\n\n\n\n" << proof << "\n\n\n\n";

  std::string proofString = proof.dump();
  std::ofstream proofFile("data/proof.json");
  if (proofFile.is_open()) {
      proofFile << proofString;
      proofFile.close();
      std::cout << "JSON data has been written to proof.json\n";
  } else {
      std::cerr << "Error opening file for writing proof.json\n";
  }
}
