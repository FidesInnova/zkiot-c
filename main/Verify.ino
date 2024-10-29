#include "FidesInnova.h"

void FidesInnova::Verify(int64_t g, int64_t mod) {
  int64_t n_i = 1;
  JsonArray array;
  bool verify = false;
  int64_t beta3 = 5;  // Must be a random number
  int64_t z_random = 2;

  int64_t alpha = 10;  // a random number  based on s_x
  int64_t beta1 = 22;  // int64_t beta1 = hashAndExtractLower4Bytes(s_x[8], mod);
  int64_t beta2 = 80;  // int64_t beta2 = hashAndExtractLower4Bytes(s_x[9], mod);
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
  g_m = ((mod - 1) * Polynomial::modInverse(m, mod)) % mod;
  y = Polynomial::power(g, g_m, mod);
  for (int64_t i = 1; i < m; i++) {
    K.push_back(Polynomial::power(y, i, mod));
  }

  vector<int64_t> H;
  int64_t w, g_n;
  H.push_back(1);
  g_n = ((mod - 1) / n) % mod;
  w = Polynomial::power(g, g_n, mod);
  for (int64_t i = 1; i < n; i++) {
    H.push_back(Polynomial::power(w, i, mod));
  }

  vector<int64_t> vH_x(n + 1, 0);
  vH_x[0] = (-1) % mod;
  if (vH_x[0] < 0) {
    vH_x[0] += mod;
  }
  vH_x[n] = 1;
  Polynomial::printPolynomial(vH_x, "vH(x)");

  vector<int64_t> vK_x = Polynomial::createLinearPolynomial(K[0]);
  for (size_t i = 1; i < K.size(); i++) {
    vector<int64_t> nextPoly = Polynomial::createLinearPolynomial(K[i]);
    vK_x = Polynomial::multiplyPolynomials(vK_x, nextPoly, mod);
  }
  Polynomial::printPolynomial(vK_x, "vK(x)");

  int64_t vH_beta1 = Polynomial::evaluatePolynomial(vH_x, beta1, mod);
  Serial.print("vH(beta1) = ");
  Serial.println(vH_beta1);

  int64_t vH_beta2 = Polynomial::evaluatePolynomial(vH_x, beta2, mod);
  Serial.print("vH(beta2) = ");
  Serial.println(vH_beta2);

  vector<int64_t> poly_beta1 = { beta1 };
  vector<int64_t> poly_beta2 = { beta2 };

  vector<int64_t> poly_pi_a = Polynomial::multiplyPolynomials(Polynomial::subtractPolynomials(rowA_x, poly_beta2, mod), Polynomial::subtractPolynomials(colA_x, poly_beta1, mod), mod);
  vector<int64_t> poly_pi_b = Polynomial::multiplyPolynomials(Polynomial::subtractPolynomials(rowB_x, poly_beta2, mod), Polynomial::subtractPolynomials(colB_x, poly_beta1, mod), mod);
  vector<int64_t> poly_pi_c = Polynomial::multiplyPolynomials(Polynomial::subtractPolynomials(rowC_x, poly_beta2, mod), Polynomial::subtractPolynomials(colC_x, poly_beta1, mod), mod);
  Polynomial::printPolynomial(poly_pi_a, "poly_pi_a(x)");
  Polynomial::printPolynomial(poly_pi_b, "poly_pi_b(x)");
  Polynomial::printPolynomial(poly_pi_c, "poly_pi_c(x)");

  int64_t poly_etaA_vH_B2_vH_B1 = (etaA * vH_beta2 * vH_beta1) % mod;
  int64_t poly_etaB_vH_B2_vH_B1 = (etaB * vH_beta2 * vH_beta1) % mod;
  int64_t poly_etaC_vH_B2_vH_B1 = (etaC * vH_beta2 * vH_beta1) % mod;

  vector<int64_t> poly_sig_a = Polynomial::multiplyPolynomialByNumber(valA_x, poly_etaA_vH_B2_vH_B1, mod);
  vector<int64_t> poly_sig_b = Polynomial::multiplyPolynomialByNumber(valB_x, poly_etaB_vH_B2_vH_B1, mod);
  vector<int64_t> poly_sig_c = Polynomial::multiplyPolynomialByNumber(valC_x, poly_etaC_vH_B2_vH_B1, mod);
  Polynomial::printPolynomial(poly_sig_a, "poly_sig_a");
  Polynomial::printPolynomial(poly_sig_b, "poly_sig_b");
  Polynomial::printPolynomial(poly_sig_c, "poly_sig_c");

  vector<int64_t> a_x = Polynomial::addPolynomials(Polynomial::addPolynomials(Polynomial::multiplyPolynomials(poly_sig_a, Polynomial::multiplyPolynomials(poly_pi_b, poly_pi_c, mod), mod), Polynomial::multiplyPolynomials(poly_sig_b, Polynomial::multiplyPolynomials(poly_pi_a, poly_pi_c, mod), mod), mod), Polynomial::multiplyPolynomials(poly_sig_c, Polynomial::multiplyPolynomials(poly_pi_a, poly_pi_b, mod), mod), mod);

  vector<int64_t> b_x = Polynomial::multiplyPolynomials(Polynomial::multiplyPolynomials(poly_pi_a, poly_pi_b, mod), poly_pi_c, mod);
  vector<int64_t> r_alpha_x = Polynomial::calculatePolynomial_r_alpha_x(alpha, n, mod);
  vector<int64_t> etaA_z_hatA_x = Polynomial::multiplyPolynomialByNumber(z_hatA, etaA, mod);
  vector<int64_t> etaB_z_hatB_x = Polynomial::multiplyPolynomialByNumber(z_hatB, etaB, mod);
  vector<int64_t> etaC_z_hatC_x = Polynomial::multiplyPolynomialByNumber(z_hatC, etaC, mod);
  vector<int64_t> Sum_M_eta_M_z_hat_M_x = Polynomial::addPolynomials(Polynomial::addPolynomials(etaA_z_hatA_x, etaB_z_hatB_x, mod), etaC_z_hatC_x, mod);

  int64_t t = n_i + 1;
  // // Calculate the known value from n^2 - n
  // int64_t n_val = (n * n - n) / 2;
  // // Calculate the value of (t^2 - t) / 2 mod
  // int target = n_val - m;
  // // Iterate through possible values of t (0 to mod-1)
  // for (int64_t i = 0; i < mod; i++) {
  //   // Calculate (t^2 - t) / 2 in mod
  //   int64_t t_val = (i * i - i) / 2;
  //   // Apply modulo to ensure positive results
  //   if (t_val < 0) t_val += mod;
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

  vector<int64_t> polyX_HAT_H = Polynomial::setupLagrangePolynomial(zero_to_t_for_H, zero_to_t_for_z, mod, "x_hat(h)");

  vector<int64_t> r_Sum_x = Polynomial::multiplyPolynomials(r_alpha_x, Sum_M_eta_M_z_hat_M_x, mod);
  vector<int64_t> v_H = Polynomial::expandPolynomials(zero_to_t_for_H, mod);
  vector<int64_t> z_hat_x = Polynomial::addPolynomials(Polynomial::multiplyPolynomials(w_hat_x, v_H, mod), polyX_HAT_H, mod);

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
  // int64_t ComP_AHP_x = (Polynomial::power(Com1_AHP_x, eta_w_hat, mod) * (Polynomial::power(Com2_AHP_x, eta_z_hatA, mod) * (Polynomial::power(Com3_AHP_x, eta_z_hatB, mod) * (Polynomial::power(Com4_AHP_x, eta_z_hatC, mod) * (Polynomial::power(Com5_AHP_x, eta_h_0_x, mod) * (Polynomial::power(Com6_AHP_x, eta_s_x, mod) * (Polynomial::power(Com7_AHP_x, eta_g_1_x, mod) * (Polynomial::power(Com8_AHP_x, eta_h_1_x, mod) * (Polynomial::power(Com9_AHP_x, eta_g_2_x, mod) * (Polynomial::power(Com10_AHP_x, eta_h_2_x, mod) * (Polynomial::power(Com11_AHP_x, eta_g_3_x, mod) * Polynomial::power(Com12_AHP_x, eta_h_3_x, mod)) % mod) % mod) % mod) % mod) % mod) % mod) % mod) % mod) % mod) % mod);

  int64_t ComP_AHP_x = Polynomial::KZG_Commitment(ck, p_x, mod);

  Polynomial::printPolynomial(a_x, "a_x");
  Polynomial::printPolynomial(b_x, "b_x");
  Polynomial::printPolynomial(g_3_x, "g_3_x");
  Serial.print("beta3 = ");
  Serial.println(beta3);
  Serial.print("sigma3 = ");
  Serial.println(sigma3);

  int64_t eq11 = (Polynomial::evaluatePolynomial(h_3_x, beta3, mod) * Polynomial::evaluatePolynomial(vK_x, beta3, mod)) % mod;
  int64_t eq12 = (Polynomial::evaluatePolynomial(a_x, beta3, mod) - ((Polynomial::evaluatePolynomial(b_x, beta3, mod) * (beta3 * Polynomial::evaluatePolynomial(g_3_x, beta3, mod) + (sigma3 * Polynomial::modInverse(m, mod)) % mod)))) % mod;
  eq12 %= mod;
  if (eq12 < 0) eq12 += mod;
  Serial.println("");
  Serial.print("eq1: ");
  Serial.print(eq11);
  Serial.print(" = ");
  Serial.println(eq12);

  int64_t eq21 = (Polynomial::evaluatePolynomial(r_alpha_x, beta2, mod) * sigma3) % mod;
  int64_t eq22 = ((Polynomial::evaluatePolynomial(h_2_x, beta2, mod) * Polynomial::evaluatePolynomial(vH_x, beta2, mod)) % mod + (beta2 * Polynomial::evaluatePolynomial(g_2_x, beta2, mod)) % mod + (sigma2 * Polynomial::modInverse(n, mod)) % mod) % mod;
  eq12 %= mod;
  if (eq12 < 0) eq12 += mod;
  Serial.print("eq2: ");
  Serial.print(eq21);
  Serial.print(" = ");
  Serial.println(eq22);

  int64_t eq31 = (Polynomial::evaluatePolynomial(s_x, beta1, mod) + Polynomial::evaluatePolynomial(r_alpha_x, beta1, mod) * Polynomial::evaluatePolynomial(Sum_M_eta_M_z_hat_M_x, beta1, mod) - sigma2 * Polynomial::evaluatePolynomial(z_hat_x, beta1, mod));
  int64_t eq32 = (Polynomial::evaluatePolynomial(h_1_x, beta1, mod) * Polynomial::evaluatePolynomial(vH_x, beta1, mod) + beta1 * Polynomial::evaluatePolynomial(g_1_x, beta1, mod) + sigma1 * Polynomial::modInverse(n, mod)) % mod;
  eq31 %= mod;
  if (eq31 < 0) eq31 += mod;
  Serial.print("eq3: ");
  Serial.print(eq31);
  Serial.print(" = ");
  Serial.println(eq32);

  int64_t eq41 = (Polynomial::evaluatePolynomial(z_hatA, beta1, mod) * Polynomial::evaluatePolynomial(z_hatB, beta1, mod) - Polynomial::evaluatePolynomial(z_hatC, beta1, mod));
  int64_t eq42 = (Polynomial::evaluatePolynomial(h_0_x, beta1, mod) * Polynomial::evaluatePolynomial(vH_x, beta1, mod)) % mod;
  eq41 %= mod;
  if (eq41 < 0) eq41 += mod;
  Serial.print("eq4: ");
  Serial.print(eq41);
  Serial.print(" = ");
  Serial.println(eq42);

  int64_t eq51Buf = (ComP_AHP_x - (g * y_prime));
  eq51Buf %= mod;
  if (eq51Buf < 0) {
    eq51Buf += mod;
  }
  int64_t eq51 = Polynomial::e_func(eq51Buf, g, g, mod);

  int64_t eq52BufP2 = (vk - (g * z_random)) % mod;
  if (eq52BufP2 < 0) {
    eq52BufP2 += mod;
  }
  int64_t eq52 = Polynomial::e_func(p_17_AHP, eq52BufP2, g, mod);
  Serial.print("eq5: ");
  Serial.print(eq51);
  Serial.print(" = ");
  Serial.println(eq52);

  if (eq11 == eq12 && eq21 == eq22 && eq31 == eq32 && eq41 == eq42 && eq51 == eq52) {
    verify = true;
  }

  Serial.println("");
  if (verify) {
    Serial.println("verify!!!!!!!!!!");
  }
}
