/*
  Read the setup.json, commitment.json, and proof.json and use them to verify the code execution 
  

  setup.json
  {
    "Class":  32-bit Integer,
    "ck": 64-bit Integer Array,
    "vk": 64-bit Integer
  }


  commitment.json
  {
    "CommitmentID":  64-bit,
    "Class": 32-bit Integer,
    "IoT_Manufacturer_Name": String,
    "IoT_Device_Name": String,
    "Device_Hardware_Version": float,
    "Firmware_Version": float,
    "Lines": 64-bit Array,
    
    "m": 64-bit Integer,
    "n": 64-bit Integer,
    
    // PFR Commitment   
    "ComRowA": 64-bit Integer,
    "ComColA": 64-bit Integer,
    "ComValA": 64-bit Integer,
    "ComRowB": 64-bit Integer,
    "ComColB": 64-bit Integer,
    "ComValB": 64-bit Integer,
    "ComRowC": 64-bit Integer,
    "ComColC": 64-bit Integer,
    "ComValC": 64-bit Integer,
    "RowA": 64-bit Integer,
    "ColA": 64-bit Integer,
    "ValA": 64-bit Integer,
    "RowB": 64-bit Integer,
    "ColB": 64-bit Integer,
    "ValB": 64-bit Integer,
    "RowC": 64-bit Integer,
    "ColC": 64-bit Integer,
    "ValC": 64-bit Integer,

    // AHP Commitment
    "ComRow'A": 64-bit Integer,
    "ComCol'A": 64-bit Integer,
    "ComVal'A": 64-bit Integer,
    "ComRow'B": 64-bit Integer,
    "ComCol'B": 64-bit Integer,
    "ComVal'B": 64-bit Integer,
    "ComRow'C": 64-bit Integer,
    "ComCol'C": 64-bit Integer,
    "ComVal'C": 64-bit Integer,
    "Row'A": 64-bit Integer,
    "Col'A": 64-bit Integer,
    "Val'A": 64-bit Integer,
    "Row'B": 64-bit Integer,
    "Col'B": 64-bit Integer,
    "Val'B": 64-bit Integer,
    "Row'C": 64-bit Integer,
    "Col'C": 64-bit Integer,
    "Val'C": 64-bit Integer,
    
    "Curve": String,
    "PolynomialCommitment": String
  }


  proof.json
  {
    "CommitmentID": 64-bit,
    "DeviceEncodedID": String,
    "Input": String,
    "Output": String,
    "P1AHP": 64-bit Integer,
    "P2AHP": 64-bit Array,
    "P3AHP": 64-bit Array,
    "P4AHP": 64-bit Array,
    "P5AHP": 64-bit Array,
    "P6AHP": 64-bit Array,
    "P7AHP": 64-bit Array,
    "P8AHP": 64-bit Array,
    "P9AHP": 64-bit Array],
    "P10AHP": 64-bit Integer,
    "P11AHP": 64-bit Array,
    "P12AHP": 64-bit Array,
    "P13AHP": 64-bit Integer,
    "P14AHP": 64-bit Array,
    "P15AHP": 64-bit Array,
    "P16AHP": 64-bit Array,    
    "Protocol": String
  }


*/

#include "fidesinnova.h"

void fidesinnova::verifier(int64_t g, int64_t p) {
  int64_t n_i = 1;
  JsonArray array;
  bool verify = false;
  int64_t beta3 = 5;  // Must be a random number
  int64_t z_random = 2;

  int64_t alpha = 10;  // a random number  based on s_x
  int64_t beta1 = 22;  // int64_t beta1 = hashAndExtractLower4Bytes(s_x[8], p);
  int64_t beta2 = 80;  // int64_t beta2 = hashAndExtractLower4Bytes(s_x[9], p);
  int64_t etaA = 2;    // a random number  based on s_x
  int64_t etaB = 30;   // a random number  based on s_x
  int64_t etaC = 100;  // a random number  based on s_x


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

  /*********************************  Read Setup  *********************************/
  String setup = readFile("/setup.json");
  DynamicJsonDocument jsonSetup(4096);  // Create a DynamicJsonDocument with a buffer size
  DeserializationError error = deserializeJson(jsonSetup, setup);

  if (error) {
    Serial.println("Failed to parse JSON");
    return;
  }

  vector<int64_t> ck;
  array = jsonSetup["ck"];
  for (JsonVariant v : array) {
    ck.push_back(v.as<int64_t>());
  }
  int64_t vk = jsonSetup["vk"][0].as<int64_t>();
  /*********************************  Read Setup  *********************************/

  /*******************************  Read Commitment  ******************************/
  String commitment = readFile("/commitment.json");
  DynamicJsonDocument jsonCommitment(4096);  // Create a DynamicJsonDocument with a buffer size
  deserializeJson(jsonCommitment, commitment);
  int64_t m = jsonCommitment["m"][0].as<int64_t>();
  int64_t n = jsonCommitment["n"][0].as<int64_t>();
  array = jsonCommitment["RowA"];
  vector<int64_t> rowA_x;
  for (JsonVariant v : array) {
    rowA_x.push_back(v.as<int64_t>());
  }
  Polynomial::printPolynomial(rowA_x, "rowA_x");

  array = jsonCommitment["ColA"];
  vector<int64_t> colA_x;
  for (JsonVariant v : array) {
    colA_x.push_back(v.as<int64_t>());
  }
  array = jsonCommitment["ValA"];
  vector<int64_t> valA_x;
  for (JsonVariant v : array) {
    valA_x.push_back(v.as<int64_t>());
  }
  array = jsonCommitment["RowB"];
  vector<int64_t> rowB_x;
  for (JsonVariant v : array) {
    rowB_x.push_back(v.as<int64_t>());
  }

  array = jsonCommitment["ColB"];
  vector<int64_t> colB_x;
  for (JsonVariant v : array) {
    colB_x.push_back(v.as<int64_t>());
  }
  array = jsonCommitment["ValB"];
  vector<int64_t> valB_x;
  for (JsonVariant v : array) {
    valB_x.push_back(v.as<int64_t>());
  }
  array = jsonCommitment["RowC"];
  vector<int64_t> rowC_x;
  for (JsonVariant v : array) {
    rowC_x.push_back(v.as<int64_t>());
  }
  array = jsonCommitment["ColC"];
  vector<int64_t> colC_x;
  for (JsonVariant v : array) {
    colC_x.push_back(v.as<int64_t>());
  }
  array = jsonCommitment["ValC"];
  vector<int64_t> valC_x;
  for (JsonVariant v : array) {
    valC_x.push_back(v.as<int64_t>());
  }
  /*******************************  Read Commitment  ******************************/

  /*********************************  Read Proof  *********************************/
  String proof = readFile("/proof.json");
  DynamicJsonDocument jsonProof(2048);  // Create a DynamicJsonDocument with a buffer size
  deserializeJson(jsonProof, proof);

  int64_t sigma1 = jsonProof["P1AHP"][0].as<int64_t>();
  array = jsonProof["P2AHP"];
  vector<int64_t> w_hat_x;
  for (JsonVariant v : array) {
    w_hat_x.push_back(v.as<int64_t>());
  }
  array = jsonProof["P3AHP"];
  vector<int64_t> z_hatA;
  for (JsonVariant v : array) {
    z_hatA.push_back(v.as<int64_t>());
  }
  array = jsonProof["P4AHP"];
  vector<int64_t> z_hatB;
  for (JsonVariant v : array) {
    z_hatB.push_back(v.as<int64_t>());
  }
  array = jsonProof["P5AHP"];
  vector<int64_t> z_hatC;
  for (JsonVariant v : array) {
    z_hatC.push_back(v.as<int64_t>());
  }
  array = jsonProof["P6AHP"];
  vector<int64_t> h_0_x;
  for (JsonVariant v : array) {
    h_0_x.push_back(v.as<int64_t>());
  }
  array = jsonProof["P7AHP"];
  vector<int64_t> s_x;
  for (JsonVariant v : array) {
    s_x.push_back(v.as<int64_t>());
  }
  array = jsonProof["P8AHP"];
  vector<int64_t> g_1_x;
  for (JsonVariant v : array) {
    g_1_x.push_back(v.as<int64_t>());
  }
  array = jsonProof["P9AHP"];
  vector<int64_t> h_1_x;
  for (JsonVariant v : array) {
    h_1_x.push_back(v.as<int64_t>());
  }
  int64_t sigma2 = jsonProof["P10AHP"][0].as<int64_t>();
  array = jsonProof["P11AHP"];
  vector<int64_t> g_2_x;
  for (JsonVariant v : array) {
    g_2_x.push_back(v.as<int64_t>());
  }
  array = jsonProof["P12AHP"];
  vector<int64_t> h_2_x;
  for (JsonVariant v : array) {
    h_2_x.push_back(v.as<int64_t>());
  }
  int64_t sigma3 = jsonProof["P13AHP"][0].as<int64_t>();
  array = jsonProof["P14AHP"];
  vector<int64_t> g_3_x;
  for (JsonVariant v : array) {
    g_3_x.push_back(v.as<int64_t>());
  }
  array = jsonProof["P15AHP"];
  vector<int64_t> h_3_x;
  for (JsonVariant v : array) {
    h_3_x.push_back(v.as<int64_t>());
  }
  int64_t y_prime = jsonProof["P16AHP"][0].as<int64_t>();
  int64_t p_17_AHP = jsonProof["P17AHP"][0].as<int64_t>();

  array = jsonProof["P18AHP"];
  vector<int64_t> z;
  for (JsonVariant v : array) {
    z.push_back(v.as<int64_t>());
  }
  /*********************************  Read Proof  *********************************/
  Serial.print("n: ");
  Serial.println(n);
  Serial.print("m: ");
  Serial.println(m);
  Serial.print("g: ");
  Serial.println(g);

  int64_t y, g_m;
  vector<int64_t> K;
  K.push_back(1);
  g_m = ((p - 1) * Polynomial::pInverse(m, p)) % p;
  y = Polynomial::power(g, g_m, p);
  for (int64_t i = 1; i < m; i++) {
    K.push_back(Polynomial::power(y, i, p));
  }

  vector<int64_t> H;
  int64_t w, g_n;
  H.push_back(1);
  g_n = ((p - 1) / n) % p;
  w = Polynomial::power(g, g_n, p);
  for (int64_t i = 1; i < n; i++) {
    H.push_back(Polynomial::power(w, i, p));
  }

  vector<int64_t> vH_x(n + 1, 0);
  vH_x[0] = (-1) % p;
  if (vH_x[0] < 0) {
    vH_x[0] += p;
  }
  vH_x[n] = 1;
  Polynomial::printPolynomial(vH_x, "vH(x)");

  vector<int64_t> vK_x = Polynomial::createLinearPolynomial(K[0]);
  for (size_t i = 1; i < K.size(); i++) {
    vector<int64_t> nextPoly = Polynomial::createLinearPolynomial(K[i]);
    vK_x = Polynomial::multiplyPolynomials(vK_x, nextPoly, p);
  }
  Polynomial::printPolynomial(vK_x, "vK(x)");

  int64_t vH_beta1 = Polynomial::evaluatePolynomial(vH_x, beta1, p);
  Serial.print("vH(beta1) = ");
  Serial.println(vH_beta1);

  int64_t vH_beta2 = Polynomial::evaluatePolynomial(vH_x, beta2, p);
  Serial.print("vH(beta2) = ");
  Serial.println(vH_beta2);

  vector<int64_t> poly_beta1 = { beta1 };
  vector<int64_t> poly_beta2 = { beta2 };

  vector<int64_t> poly_pi_a = Polynomial::multiplyPolynomials(Polynomial::subtractPolynomials(rowA_x, poly_beta2, p), Polynomial::subtractPolynomials(colA_x, poly_beta1, p), p);
  vector<int64_t> poly_pi_b = Polynomial::multiplyPolynomials(Polynomial::subtractPolynomials(rowB_x, poly_beta2, p), Polynomial::subtractPolynomials(colB_x, poly_beta1, p), p);
  vector<int64_t> poly_pi_c = Polynomial::multiplyPolynomials(Polynomial::subtractPolynomials(rowC_x, poly_beta2, p), Polynomial::subtractPolynomials(colC_x, poly_beta1, p), p);
  Polynomial::printPolynomial(poly_pi_a, "poly_pi_a(x)");
  Polynomial::printPolynomial(poly_pi_b, "poly_pi_b(x)");
  Polynomial::printPolynomial(poly_pi_c, "poly_pi_c(x)");

  vector<int64_t> poly_etaA_vH_B2_vH_B1 = { (etaA * vH_beta2 * vH_beta1) % p };
  vector<int64_t> poly_etaB_vH_B2_vH_B1 = { (etaB * vH_beta2 * vH_beta1) % p };
  vector<int64_t> poly_etaC_vH_B2_vH_B1 = { (etaC * vH_beta2 * vH_beta1) % p };

  vector<int64_t> poly_sig_a = Polynomial::multiplyPolynomials(poly_etaA_vH_B2_vH_B1, valA_x, p);
  vector<int64_t> poly_sig_b = Polynomial::multiplyPolynomials(poly_etaB_vH_B2_vH_B1, valB_x, p);
  vector<int64_t> poly_sig_c = Polynomial::multiplyPolynomials(poly_etaC_vH_B2_vH_B1, valC_x, p);
  Polynomial::printPolynomial(poly_sig_a, "poly_sig_a(x)");
  Polynomial::printPolynomial(poly_sig_b, "poly_sig_b(x)");
  Polynomial::printPolynomial(poly_sig_c, "poly_sig_c(x)");

  vector<int64_t> a_x = Polynomial::addPolynomials(Polynomial::addPolynomials(Polynomial::multiplyPolynomials(poly_sig_a, Polynomial::multiplyPolynomials(poly_pi_b, poly_pi_c, p), p), Polynomial::multiplyPolynomials(poly_sig_b, Polynomial::multiplyPolynomials(poly_pi_a, poly_pi_c, p), p), p), Polynomial::multiplyPolynomials(poly_sig_c, Polynomial::multiplyPolynomials(poly_pi_a, poly_pi_b, p), p), p);

  vector<int64_t> b_x = Polynomial::multiplyPolynomials(Polynomial::multiplyPolynomials(poly_pi_a, poly_pi_b, p), poly_pi_c, p);
  vector<int64_t> r_alpha_x = Polynomial::calculatePolynomial_r_alpha_x(alpha, n, p);
  vector<int64_t> etaA_z_hatA_x = Polynomial::multiplyPolynomialByNumber(z_hatA, etaA, p);
  vector<int64_t> etaB_z_hatB_x = Polynomial::multiplyPolynomialByNumber(z_hatB, etaB, p);
  vector<int64_t> etaC_z_hatC_x = Polynomial::multiplyPolynomialByNumber(z_hatC, etaC, p);
  vector<int64_t> Sum_M_eta_M_z_hat_M_x = Polynomial::addPolynomials(Polynomial::addPolynomials(etaA_z_hatA_x, etaB_z_hatB_x, p), etaC_z_hatC_x, p);

  int64_t t = n_i + 1;
  // // Calculate the known value from n^2 - n
  // int64_t n_val = (n * n - n) / 2;
  // // Calculate the value of (t^2 - t) / 2 p
  // int target = n_val - m;
  // // Iterate through possible values of t (0 to p-1)
  // for (int64_t i = 0; i < p; i++) {
  //   // Calculate (t^2 - t) / 2 in p
  //   int64_t t_val = (i * i - i) / 2;
  //   // Apply pulo to ensure positive results
  //   if (t_val < 0) t_val += p;
  //   if (t_val == target) {
  //     t = i;
  //   }
  // }
  Serial.print("t = ");
  Serial.println(t);

  vector<int64_t> zero_to_t_for_H;
  vector<int64_t> zero_to_t_for_z;
  for (int64_t i = 0; i < t; i++) {
    zero_to_t_for_H.push_back(H[i]);
    zero_to_t_for_z.push_back(z[i]);
  }

  vector<int64_t> polyX_HAT_H = Polynomial::setupLagrangePolynomial(zero_to_t_for_H, zero_to_t_for_z, p, "x_hat(h)");

  vector<int64_t> r_Sum_x = Polynomial::multiplyPolynomials(r_alpha_x, Sum_M_eta_M_z_hat_M_x, p);
  vector<int64_t> v_H = Polynomial::expandPolynomials(zero_to_t_for_H, p);
  vector<int64_t> z_hat_x = Polynomial::addPolynomials(Polynomial::multiplyPolynomials(w_hat_x, v_H, p), polyX_HAT_H, p);

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
  // int64_t ComP_AHP_x = (Polynomial::power(Com1_AHP_x, eta_w_hat, p) * (Polynomial::power(Com2_AHP_x, eta_z_hatA, p) * (Polynomial::power(Com3_AHP_x, eta_z_hatB, p) * (Polynomial::power(Com4_AHP_x, eta_z_hatC, p) * (Polynomial::power(Com5_AHP_x, eta_h_0_x, p) * (Polynomial::power(Com6_AHP_x, eta_s_x, p) * (Polynomial::power(Com7_AHP_x, eta_g_1_x, p) * (Polynomial::power(Com8_AHP_x, eta_h_1_x, p) * (Polynomial::power(Com9_AHP_x, eta_g_2_x, p) * (Polynomial::power(Com10_AHP_x, eta_h_2_x, p) * (Polynomial::power(Com11_AHP_x, eta_g_3_x, p) * Polynomial::power(Com12_AHP_x, eta_h_3_x, p)) % p) % p) % p) % p) % p) % p) % p) % p) % p) % p);

  int64_t ComP_AHP_x = Polynomial::KZG_Commitment(ck, p_x, p);

  Polynomial::printPolynomial(a_x, "a_x");
  Polynomial::printPolynomial(b_x, "b_x");
  Polynomial::printPolynomial(g_3_x, "g_3_x");
  Serial.print("beta3 = ");
  Serial.println(beta3);
  Serial.print("sigma3 = ");
  Serial.println(sigma3);

  int64_t eq11 = (Polynomial::evaluatePolynomial(h_3_x, beta3, p) * Polynomial::evaluatePolynomial(vK_x, beta3, p)) % p;
  int64_t eq12 = (Polynomial::evaluatePolynomial(a_x, beta3, p) - ((Polynomial::evaluatePolynomial(b_x, beta3, p) * (beta3 * Polynomial::evaluatePolynomial(g_3_x, beta3, p) + (sigma3 * Polynomial::pInverse(m, p)) % p)))) % p;
  eq12 %= p;
  if (eq12 < 0) eq12 += p;
  Serial.print(eq11);
  Serial.print(" = ");
  Serial.println(eq12);

  int64_t eq21 = (Polynomial::evaluatePolynomial(r_alpha_x, beta2, p) * sigma3) % p;
  int64_t eq22 = ((Polynomial::evaluatePolynomial(h_2_x, beta2, p) * Polynomial::evaluatePolynomial(vH_x, beta2, p)) % p + (beta2 * Polynomial::evaluatePolynomial(g_2_x, beta2, p)) % p + (sigma2 * Polynomial::pInverse(n, p)) % p) % p;
  eq12 %= p;
  if (eq12 < 0) eq12 += p;
  Serial.print(eq21);
  Serial.print(" = ");
  Serial.println(eq22);

  int64_t eq31 = (Polynomial::evaluatePolynomial(s_x, beta1, p) + Polynomial::evaluatePolynomial(r_alpha_x, beta1, p) * Polynomial::evaluatePolynomial(Sum_M_eta_M_z_hat_M_x, beta1, p) - sigma2 * Polynomial::evaluatePolynomial(z_hat_x, beta1, p));
  int64_t eq32 = (Polynomial::evaluatePolynomial(h_1_x, beta1, p) * Polynomial::evaluatePolynomial(vH_x, beta1, p) + beta1 * Polynomial::evaluatePolynomial(g_1_x, beta1, p) + sigma1 * Polynomial::pInverse(n, p)) % p;
  eq31 %= p;
  if (eq31 < 0) eq31 += p;
  Serial.print(eq31);
  Serial.print(" = ");
  Serial.println(eq32);

  int64_t eq41 = (Polynomial::evaluatePolynomial(z_hatA, beta1, p) * Polynomial::evaluatePolynomial(z_hatB, beta1, p) - Polynomial::evaluatePolynomial(z_hatC, beta1, p));
  int64_t eq42 = (Polynomial::evaluatePolynomial(h_0_x, beta1, p) * Polynomial::evaluatePolynomial(vH_x, beta1, p)) % p;
  eq41 %= p;
  if (eq41 < 0) eq41 += p;
  Serial.print(eq41);
  Serial.print(" = ");
  Serial.println(eq42);

  int64_t eq51Buf = (ComP_AHP_x - (g * y_prime));
  eq51Buf %= p;
  if (eq51Buf < 0) {
    eq51Buf += p;
  }
  int64_t eq51 = Polynomial::e_func(eq51Buf, g, g, p);

  int64_t eq52BufP2 = (vk - (g * z_random)) % p;
  if (eq52BufP2 < 0) {
    eq52BufP2 += p;
  }
  int64_t eq52 = Polynomial::e_func(p_17_AHP, eq52BufP2, g, p);
  Serial.print("eq51 = ");
  Serial.println(eq51);
  Serial.print("eq52 = ");
  Serial.println(eq52);

  if (eq11 == eq12 && eq21 == eq22 && eq31 == eq32 && eq41 == eq42 && eq51 == eq52) {
    verify = true;
  }

  Serial.println("");
  if (verify) {
    Serial.println("verify!!!!!!!!!!");
  }
}
