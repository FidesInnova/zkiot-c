#include "FidesInnova.h"

void FidesInnova::Verify(int64_t g, int64_t mod) {
  JsonArray array;
  bool verify = false;
  int64_t beta3 = 5;  // Must be a random number

  /*********************************  Read Setup  *********************************/
  String setup = readFile("/setup.json");
  DynamicJsonDocument jsonSetup(2048);  // Create a DynamicJsonDocument with a buffer size
  deserializeJson(jsonSetup, setup);
  vector<int64_t> ck;
  array = jsonSetup["ck"];
  for (JsonVariant v : array) {
    ck.push_back(v.as<int64_t>());
  }
  int64_t vk = jsonSetup["vk"];
  /*********************************  Read Setup  *********************************/

  /*******************************  Read Commitment  ******************************/
  String commitment = readFile("/commitment.json");
  DynamicJsonDocument jsonCommitment(2048);  // Create a DynamicJsonDocument with a buffer size
  deserializeJson(jsonCommitment, commitment);
  int64_t m = jsonCommitment["m"];
  int64_t n = jsonCommitment["n"];
  array = jsonSetup["RowA"];
  vector<int64_t> rowA_x;
  for (JsonVariant v : array) {
    rowA_x.push_back(v.as<int64_t>());
  }
  array = jsonSetup["ColA"];
  vector<int64_t> colA_x;
  for (JsonVariant v : array) {
    colA_x.push_back(v.as<int64_t>());
  }
  array = jsonSetup["ValA"];
  vector<int64_t> valA_x;
  for (JsonVariant v : array) {
    valA_x.push_back(v.as<int64_t>());
  }
  array = jsonSetup["RowB"];
  vector<int64_t> rowB_x;
  for (JsonVariant v : array) {
    rowB_x.push_back(v.as<int64_t>());
  }
  array = jsonSetup["ColB"];
  vector<int64_t> colB_x;
  for (JsonVariant v : array) {
    colB_x.push_back(v.as<int64_t>());
  }
  array = jsonSetup["ValB"];
  vector<int64_t> valB_x;
  for (JsonVariant v : array) {
    valB_x.push_back(v.as<int64_t>());
  }
  array = jsonSetup["RowC"];
  vector<int64_t> rowC_x;
  for (JsonVariant v : array) {
    rowC_x.push_back(v.as<int64_t>());
  }
  array = jsonSetup["ColC"];
  vector<int64_t> colC_x;
  for (JsonVariant v : array) {
    colC_x.push_back(v.as<int64_t>());
  }
  array = jsonSetup["ValC"];
  vector<int64_t> valC_x;
  for (JsonVariant v : array) {
    valC_x.push_back(v.as<int64_t>());
  }
  /*******************************  Read Commitment  ******************************/

  /*********************************  Read Proof  *********************************/
  String proof = readFile("/proof.json");
  DynamicJsonDocument jsonProof(2048);  // Create a DynamicJsonDocument with a buffer size
  deserializeJson(jsonProof, proof);

  int64_t sigma1 = jsonProof["P1AHP"];
  array = jsonSetup["P2AHP"];
  vector<int64_t> w_hat_x;
  for (JsonVariant v : array) {
    w_hat_x.push_back(v.as<int64_t>());
  }
  array = jsonSetup["P3AHP"];
  vector<int64_t> z_hatA;
  for (JsonVariant v : array) {
    z_hatA.push_back(v.as<int64_t>());
  }
  array = jsonSetup["P4AHP"];
  vector<int64_t> z_hatB;
  for (JsonVariant v : array) {
    z_hatB.push_back(v.as<int64_t>());
  }
  array = jsonSetup["P5AHP"];
  vector<int64_t> z_hatC;
  for (JsonVariant v : array) {
    z_hatC.push_back(v.as<int64_t>());
  }
  array = jsonSetup["P6AHP"];
  vector<int64_t> h_0_x;
  for (JsonVariant v : array) {
    h_0_x.push_back(v.as<int64_t>());
  }
  array = jsonSetup["P7AHP"];
  vector<int64_t> s_x;
  for (JsonVariant v : array) {
    s_x.push_back(v.as<int64_t>());
  }
  array = jsonSetup["P8AHP"];
  vector<int64_t> g_1_x;
  for (JsonVariant v : array) {
    g_1_x.push_back(v.as<int64_t>());
  }
  array = jsonSetup["P9AHP"];
  vector<int64_t> h_1_x;
  for (JsonVariant v : array) {
    h_1_x.push_back(v.as<int64_t>());
  }
  int64_t sigma2 = jsonProof["P10AHP"];
  
  array = jsonSetup["P11AHP"];
  vector<int64_t> g_2_x;
  for (JsonVariant v : array) {
    g_2_x.push_back(v.as<int64_t>());
  }
  array = jsonSetup["P12AHP"];
  vector<int64_t> h_2_x;
  for (JsonVariant v : array) {
    h_2_x.push_back(v.as<int64_t>());
  }
  int64_t sigma3 = jsonProof["P13AHP"];
  array = jsonSetup["P14AHP"];
  vector<int64_t> g_3_x;
  for (JsonVariant v : array) {
    g_3_x.push_back(v.as<int64_t>());
  }
  array = jsonSetup["P15AHP"];
  vector<int64_t> h_3_x;
  for (JsonVariant v : array) {
    h_3_x.push_back(v.as<int64_t>());
  }
  int64_t y_prime = jsonProof["P16AHP"];
  array = jsonSetup["P17AHP"];
  vector<int64_t> p_17_AHP;
  for (JsonVariant v : array) {
    p_17_AHP.push_back(v.as<int64_t>());
  }
  /*********************************  Read Proof  *********************************/


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
  Serial.print(eq41);
  Serial.print(" = ");
  Serial.println(eq42);

  int64_t eq51Buf = (ComP_AHP_x * Polynomial::modInverse(Polynomial::power(g, y_prime, mod), mod)) % mod;
  Serial.print("ef1 = ");
  Serial.println(eq51Buf);
  int64_t eq51 = Polynomial::e_func(eq51Buf, g, g, mod);
  int64_t eq52BufP2 = (vk * Polynomial::modInverse(Polynomial::power(g, z_random, mod), mod)) % mod;
  int64_t eq52 = Polynomial::e_func(p_17_AHP, eq52BufP2, g, mod);
  // vector<int64_t> eq52BufP1;
  // eq52BufP1.push_back(p_17_AHP);
  // int64_t eq52 = (Polynomial::LagrangePolynomial(eq52BufP1, g, mod) * Polynomial::LagrangePolynomial(eq52BufP2, g, mod)) % mod;
  Serial.print("eq51 = ");
  Serial.println(eq51);
  Serial.print("eq52 = ");
  Serial.println(eq52);

  if (eq11 == eq12 && eq21 == eq22 && eq31 == eq32 && eq41 == eq42) {
    verify = true;
  }

  Serial.println("");
  if (verify) {
    Serial.println("verify!!!!!!!!!!");
  }
}
