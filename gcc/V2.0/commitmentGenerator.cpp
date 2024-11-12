/*
  Read the directory that containing program.s, device_config.json, and setup.json

  program.s
    .file "program.cpp"
    .option nopic
    .attribute arch, "rv32i2p1_m2p0_a2p1_f2p2_d2p2_c2p0_zicsr2p0_zifencei2p0"
    .attribute unaligned_access, 0
    .attribute stack_align, 16
    .text
    .align  1
    .globl  main
    .type main, @function
  main:
  .LFB0:
    .cfi_startproc
    addi  sp,sp,-32
    .cfi_def_cfa_offset 32
    sw  ra,28(sp)
    sw  s0,24(sp)
    .cfi_offset 1, -4
    .cfi_offset 8, -8
    addi  s0,sp,32
    .cfi_def_cfa 8, 0
    li  a5,297
    sw  a5,-20(s0)
    j .L2
  .L3:
    lw  a5,-20(s0)
    addi  a5,a5,383
    sw  a5,-20(s0)
  .L2:
    lw  a4,-20(s0)
    li  a5,8192
    addi  a5,a5,1807
    ble a4,a5,.L3
    li  a5,0
    mv  a0,a5
    lw  ra,28(sp)
    .cfi_restore 1
    lw  s0,24(sp)
    .cfi_restore 8
    .cfi_def_cfa 2, 32
    addi  sp,sp,32
    .cfi_def_cfa_offset 0
    jr  ra
    .cfi_endproc
  .LFE0:
    .size main, .-main
    .ident  "GCC: (g04696df09) 14.2.0"
    .section  .note.GNU-stack,"",@progbits


Read the Class and code_block from device_config.json to be used in the code and pass all parameters to the program_commitment.json.
  device_config.json
  {
    "Class": 32-bit Integer,
    "IoT_Manufacturer_Name": String,
    "IoT_Device_Name": String,
    "Device_Hardware_Version": float,
    "Firmware_Version": float,
    "code_block": 64-bit Array
  }


  setup.json
  {
    "Class":  32-bit Integer,
    "ck": 64-bit Integer Array,
    "vk": 64-bit Integer
  }

*/

// #include "fidesinnova.h"

#include "polynomial.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "json.hpp"
using ordered_json = nlohmann::ordered_json;
#include <regex>
#include <sstream>
#include <unordered_map>

using namespace std;

// Map of RISC-V register aliases to x registers
std::unordered_map<std::string, int> registerMap = {
    {"zero", 0},  {"ra", 1},   {"sp", 2},   {"gp", 3},   {"tp", 4},
    {"t0", 5},    {"t1", 6},   {"t2", 7},   {"s0", 8},   {"s1", 9},
    {"a0", 10},   {"a1", 11},  {"a2", 12},  {"a3", 13},  {"a4", 14},
    {"a5", 15},   {"a6", 16},  {"a7", 17},  {"s2", 18},  {"s3", 19},
    {"s4", 20},   {"s5", 21},  {"s6", 22},  {"s7", 23},  {"s8", 24},
    {"s9", 25},   {"s10", 26}, {"s11", 27}, {"t3", 28},  {"t4", 29},
    {"t5", 30},   {"t6", 31}
};

std::vector<std::string> instructions;
int64_t Class;
string commitmentID;
string IoT_Manufacturer_Name;
string IoT_Device_Name;
float Device_Hardware_Version;
float Firmware_Version;

// Function to read JSON config file and parse lines to read from assembly file
std::pair<int64_t, int64_t> parseDeviceConfig(const std::string &configFile, nlohmann::json &config) {
  std::ifstream configFileStream(configFile, std::ifstream::binary);
  if (!configFileStream.is_open()) {
      std::cerr << "Error opening config file: " << configFile << std::endl;
      exit(EXIT_FAILURE);
  }

  configFileStream >> config;
  configFileStream.close();

  std::vector<int64_t> linesToRead;

  int64_t startLine = config["code_block"][0].get<int64_t>();
  int64_t endLine = config["code_block"][1].get<int64_t>();
  Class = config["Class"].get<int64_t>();
  commitmentID = config["commitmentID"].get<string>();
  IoT_Manufacturer_Name = config["IoT_Manufacturer_Name"].get<string>();
  IoT_Device_Name = config["IoT_Device_Name"].get<string>();
  Device_Hardware_Version = config["Device_Hardware_Version"].get<float>();
  Firmware_Version = config["Firmware_Version"].get<float>();

  return {startLine, endLine};
}

// Function to read specified lines from assembly file
std::vector<std::string> readAssemblyLines(const std::string &assemblyFile, int64_t startLine, int64_t endLine) {
  std::ifstream assemblyFileStream(assemblyFile);
  if (!assemblyFileStream.is_open()) {
      std::cerr << "Error opening assembly file: " << assemblyFile << std::endl;
      exit(EXIT_FAILURE);
  }

  std::vector<std::string> selectedLines;
  std::string line;
  int64_t currentLineNumber = 1;

  while (std::getline(assemblyFileStream, line)) {
      if (currentLineNumber >= startLine && currentLineNumber <= endLine) {
          selectedLines.push_back(line);
      }
      ++currentLineNumber;
  }

  assemblyFileStream.close();
  return selectedLines;
}

// Function to modify assembly file content and save to new file
void modifyAndSaveAssembly(const std::string &assemblyFile, const std::string &newAssemblyFile, int64_t startLine, int64_t endLine) {
  std::ifstream assemblyFileStream(assemblyFile);
  if (!assemblyFileStream.is_open()) {
      std::cerr << "Error opening assembly file: " << assemblyFile << std::endl;
      exit(EXIT_FAILURE);
  }

  std::ofstream newAssemblyFileStream(newAssemblyFile);
  if (!newAssemblyFileStream.is_open()) {
      std::cerr << "Error creating new assembly file: " << newAssemblyFile << std::endl;
      exit(EXIT_FAILURE);
  }

  std::string line;
  int64_t currentLineNumber = 1;
  int64_t index = 0;

  vector<int64_t> spaceSize(32, 4);
  vector<int64_t> rdList;
  while (std::getline(assemblyFileStream, line)) {
    // Insert variables before the specified lines
    if (currentLineNumber == startLine) {
      newAssemblyFileStream << "jal store_register_instances\n";
      newAssemblyFileStream << line << std::endl;
      std::stringstream ss(line);
      std::string opcode, rd, leftStr, rightStr;
      ss >> opcode >> rd >> leftStr >> rightStr;
      rd = Polynomial::trim(rd);
      rd = Polynomial::removeCommas(rd);
      instructions.push_back(line);
      rdList.push_back(registerMap[rd]);
      newAssemblyFileStream << "la t0, x" << std::to_string(registerMap[rd]) << "_array" << endl;
      newAssemblyFileStream << "sw x" << std::to_string(registerMap[rd]) << ", " << std::to_string(spaceSize[registerMap[rd]]) << "(t0)" << endl;
      spaceSize[registerMap[rd]] += 4;
    }
    else if(currentLineNumber > startLine && currentLineNumber <= endLine){
      newAssemblyFileStream << line << std::endl;
      std::stringstream ss(line);
      std::string opcode, rd, leftStr, rightStr;
      ss >> opcode >> rd >> leftStr >> rightStr;
      rd = Polynomial::trim(rd);
      rd = Polynomial::removeCommas(rd);
      instructions.push_back(line);
      rdList.push_back(registerMap[rd]);
      newAssemblyFileStream << "la t0, x" << std::to_string(registerMap[rd]) << "_array" << endl;
      newAssemblyFileStream << "sw x" << std::to_string(registerMap[rd]) << ", " << std::to_string(spaceSize[registerMap[rd]]) << "(t0)\n";
      spaceSize[registerMap[rd]] += 4;
    }
    else if (currentLineNumber == endLine + 1){
      for(int64_t i = 0; i < n_i; i++) {
        
      }
      newAss << rdList
      newAssemblyFileStream << line << std::endl;
    }
    else {
      newAssemblyFileStream << line << std::endl;
    }

    ++currentLineNumber;
  }

  std::string assemblyCode = "#### Subroutine Code (`store_registers.s`)\n\n";
  assemblyCode +=
  "        .data\n";
  for (int i = 0; i < 32; i++) {
    assemblyCode += "x" + std::to_string(i) + "_array:    .space " + std::to_string(spaceSize[i]) + "   # Array for x" + std::to_string(i) + "\n";
  }
  assemblyCode += "z_array:    .space " + std::to_string((n_i + n_g + 1) * 4) + "   # Array for z\n";

  assemblyCode += "\n    .text\n"
  "      .globl store_register_instances\n"
  "  store_register_instances:\n"
  "      # Store each register's value in its respective array\n""      la t0, x0_array\n"
  "      sw x0, 0(t0)            # Store x0 in x0_array\n"
  "      la t0, x1_array\n"
  "      sw x1, 0(t0)            # Store x1 in x1_array\n"
  "      la t0, x2_array\n"
  "      sw x2, 0(t0)            # Store x2 in x2_array\n"
  "      la t0, x3_array\n"
  "      sw x3, 0(t0)            # Store x3 in x3_array\n"
  "      la t0, x4_array\n"
  "      sw x4, 0(t0)            # Store x4 in x4_array\n"
  "      la t0, x5_array\n"
  "      sw x5, 0(t0)            # Store x5 in x5_array\n"
  "      la t0, x6_array\n"
  "      sw x6, 0(t0)            # Store x6 in x6_array\n"
  "      la t0, x7_array\n"
  "      sw x7, 0(t0)            # Store x7 in x7_array\n"
  "      la t0, x8_array\n"
  "      sw x8, 0(t0)            # Store x8 in x8_array\n"
  "      la t0, x9_array\n"
  "      sw x9, 0(t0)            # Store x9 in x9_array\n"
  "      la t0, x10_array\n"
  "      sw x10, 0(t0)           # Store x10 in x10_array\n"
  "      la t0, x11_array\n"
  "      sw x11, 0(t0)           # Store x11 in x11_array\n"
  "      la t0, x12_array\n"
  "      sw x12, 0(t0)           # Store x12 in x12_array\n"
  "      la t0, x13_array\n"
  "      sw x13, 0(t0)           # Store x13 in x13_array\n"
  "      la t0, x14_array\n"
  "      sw x14, 0(t0)           # Store x14 in x14_array\n"
  "      la t0, x15_array\n"
  "      sw x15, 0(t0)           # Store x15 in x15_array\n"
  "      la t0, x16_array\n"
  "      sw x16, 0(t0)           # Store x16 in x16_array\n"
  "      la t0, x17_array\n"
  "      sw x17, 0(t0)           # Store x17 in x17_array\n"
  "      la t0, x18_array\n"
  "      sw x18, 0(t0)           # Store x18 in x18_array\n"
  "      la t0, x19_array\n"
  "      sw x19, 0(t0)           # Store x19 in x19_array\n"
  "      la t0, x20_array\n"
  "      sw x20, 0(t0)           # Store x20 in x20_array\n"
  "      la t0, x21_array\n"
  "      sw x21, 0(t0)           # Store x21 in x21_array\n"
  "      la t0, x22_array\n"
  "      sw x22, 0(t0)           # Store x22 in x22_array\n"
  "      la t0, x23_array\n"
  "      sw x23, 0(t0)           # Store x23 in x23_array\n"
  "      la t0, x24_array\n"
  "      sw x24, 0(t0)           # Store x24 in x24_array\n"
  "      la t0, x25_array\n"
  "      sw x25, 0(t0)           # Store x25 in x25_array\n"
  "      la t0, x26_array\n"
  "      sw x26, 0(t0)           # Store x26 in x26_array\n"
  "      la t0, x27_array\n"
  "      sw x27, 0(t0)           # Store x27 in x27_array\n"
  "      la t0, x28_array\n"
  "      sw x28, 0(t0)           # Store x28 in x28_array\n"
  "      la t0, x29_array\n"
  "      sw x29, 0(t0)           # Store x29 in x29_array\n"
  "      la t0, x30_array\n"
  "      sw x30, 0(t0)           # Store x30 in x30_array\n"
  "      la t0, x31_array\n"
  "      sw x31, 0(t0)           # Store x31 in x31_array\n"
  "      \n"
  "      ret                            # Return from function\n";


  // Replace all instances of "{SPACE_SIZE}" with the actual value of spaceSize
  size_t pos = 0;
  while ((pos = assemblyCode.find("{SPACE_SIZE}", pos)) != std::string::npos) {
    // Convert spaceSize to a string
    // std::ostringstream oss;
    // oss << spaceSize[pos];
    std::string spaceSizeStr = std::to_string(spaceSize[pos]);

    assemblyCode.replace(pos, spaceSizeStr.length(), spaceSizeStr); // spaceSize is the length of "{SPACE_SIZE}"
    pos += spaceSizeStr.length();
  }

  newAssemblyFileStream << assemblyCode << std::endl;

  assemblyFileStream.close();
  newAssemblyFileStream.close();
}


void commitmentGenerator(const std::string &setupFile) {
  std::ifstream setupFileStream(setupFile);
  if (!setupFileStream.is_open()) {
      std::cerr << "Could not open the file!" << std::endl;
  }
  nlohmann::json setupJsonData;
  setupFileStream >> setupJsonData;
  setupFileStream.close();
  vector<int64_t> ck = setupJsonData["ck"].get<vector<int64_t>>();
  int64_t vk = setupJsonData["vk"].get<int64_t>();

  
  std::ifstream classFileStream("class.json");
  if (!classFileStream.is_open()) {
      std::cerr << "Could not open the file!" << std::endl;
  }
  nlohmann::json classJsonData;
  classFileStream >> classJsonData;
  classFileStream.close();
  int64_t n_i, n_g, m, n, p, g;
  for (const auto& item : classJsonData) {
    if (item["Class"] == Class) {
      // Number of inputs, gates, m, n, p, and g
      n_i = item["n_i"].get<int64_t>();
      n_g = item["n_g"].get<int64_t>();
      m = item["m"].get<int64_t>();
      n = item["n"].get<int64_t>();
      p = item["p"].get<int64_t>();
      g = item["g"].get<int64_t>();

      std::cout << "n_i: " << n_i << ", n_g: " << n_g << std::endl;
      break; // Stop after finding the first matching "Class"
    }
  }


 for (const auto& instr : instructions) {
    std::stringstream ss(instr);
    std::string opcode, rd, leftStr, rightStr;
    
    ss >> opcode >> rd;
    ss >> leftStr >> rightStr;
    leftStr = Polynomial::trim(leftStr);
    rightStr = Polynomial::trim(rightStr);
    leftStr = Polynomial::removeCommas(leftStr);
    rightStr = Polynomial::removeCommas(rightStr);
    // std::cout << "opcode: " << opcode << "\tleftInt: " << leftStr << "\trightInt: " << rightStr << "\n";
  }
  cout << "Number of immediate instructions (n_i): " << n_i << endl;
  cout << "Number of general instructions (n_g): " << n_g << endl;

  // Matrix order
  int64_t t;
  cout << "Matrix order: " << n << endl;

  t = n_i + 1;
  // m = (((Polynomial::power(n, 2, p) - n) / 2) - ((Polynomial::power(t, 2, p) - t) / 2)) % p;

  // Initialize matrices A, B, C
  vector<vector<int64_t>> A(n, vector<int64_t>(n, 0ll));
  vector<vector<int64_t>> B(n, vector<int64_t>(n, 0ll));
  vector<vector<int64_t>> C(n, vector<int64_t>(n, 0ll));

  // Fill matrices based on the instructions
  for (int64_t i = 0; i < n_g; i++) {
    std::stringstream ss(instructions[i]);
    std::string opcode, rd, leftStr, rightStr;
    ss >> opcode >> rd;
    int64_t R = std::stoi(rd.substr(1));

    if (opcode == "add" || opcode == "addi" || opcode == "mul") {
      ss >> leftStr >> rightStr;

      // Remove commas
      leftStr = Polynomial::removeCommas(leftStr);
      rightStr = Polynomial::removeCommas(rightStr);
      // Trim spaces
      leftStr = Polynomial::trim(leftStr);
      rightStr = Polynomial::trim(rightStr);

      int64_t leftInt, rightInt;
      
      C[1+n_i+i][1+n_i+i] = 1;

      if (opcode == "add" || opcode == "addi") {
        A[1+n_i+i][0] = 1;
        if (std::isdigit(leftStr[0])) {
          leftInt = std::stoi(leftStr);
          B[1+n_i+i][0] = leftInt;
        }
        else {
          B[1+n_i+i][1 + registerMap[leftStr]] = 1;
        }
        if(std::isdigit(rightStr[0])){
          rightInt = std::stoi(rightStr);
          B[1+n_i+i][0] = rightInt;
        }
        else {
          B[1+n_i+i][1 + registerMap[rightStr]] = 1;
        }

    } else if (opcode == "mul") {
        if (std::isdigit(leftStr[0])) {
          leftInt = std::stoi(leftStr);
          A[1+n_i+i][0] = leftInt;
        }
        else {
          A[1+n_i+i][1 + registerMap[leftStr]] = 1;
        }
        if (std::isdigit(rightStr[0])) {
          rightInt = std::stoi(rightStr);
          B[1+n_i+i][0] = rightInt;
        }
        else {
          B[1+n_i+i][1 + registerMap[rightStr]] = 1;
        }
      }
    }
    
    else {
      cout << "!!! Undefined instruction in the defiend Line range !!!\n" << opcode << endl;
      std::exit(0);
    }
  }

  Polynomial::printMatrix(A, "A");
  Polynomial::printMatrix(B, "B");
  Polynomial::printMatrix(C, "C");

  // Vector H to store powers of w
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

  // Vector K to store powers of y
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
  
  // Create a polynomial vector vH_x of size (n + 1) initialized to 0
  vector<int64_t> vH_x(n + 1, 0);
  vH_x[0] = (-1) % p;
  if (vH_x[0] < 0) {
    vH_x[0] += p;
  }
  vH_x[n] = 1;
  Polynomial::printPolynomial(vH_x, "vH(x)");

 // Create a mapping for the non-zero rows using parameters K and H
  vector<vector<int64_t>> nonZeroRowsA = Polynomial::getNonZeroRows(A);
  vector<vector<int64_t>> rowA = Polynomial::createMapping(K, H, nonZeroRowsA);
  // rowA[1].push_back(1);
  // rowA[1].push_back(135);
  // rowA[1].push_back(125);
  // rowA[1].push_back(59);
  // rowA[1].push_back(42);
  // rowA[1].push_back(1);
  Polynomial::printMapping(rowA, "row_A");
  vector<vector<int64_t>> nonZeroColsA = Polynomial::getNonZeroCols(A);
  vector<vector<int64_t>> colA = Polynomial::createMapping(K, H, nonZeroColsA);
  // colA[1].push_back(42);
  // colA[1].push_back(1);
  // colA[1].push_back(135);
  // colA[1].push_back(125);
  // colA[1].push_back(59);
  // colA[1].push_back(42);
  Polynomial::printMapping(colA, "col_A");
  vector<vector<int64_t>> valA = Polynomial::valMapping(K, H, nonZeroRowsA, nonZeroColsA, p);
  Polynomial::printMapping(valA, "val_A");

  vector<vector<int64_t>> nonZeroRowsB = Polynomial::getNonZeroRows(B);
  vector<vector<int64_t>> rowB = Polynomial::createMapping(K, H, nonZeroRowsB);
  // rowB[1].push_back(59);
  // rowB[1].push_back(1);
  // rowB[1].push_back(42);
  // rowB[1].push_back(135);
  // rowB[1].push_back(59);
  Polynomial::printMapping(rowB, "row_B");
  vector<vector<int64_t>> nonZeroColsB = Polynomial::getNonZeroCols(B);
  vector<vector<int64_t>> colB = Polynomial::createMapping(K, H, nonZeroColsB);
  // colB[1].push_back(59);
  // colB[1].push_back(42);
  // colB[1].push_back(125);
  // colB[1].push_back(1);
  // colB[1].push_back(135);
  Polynomial::printMapping(colB, "col_B");
  vector<vector<int64_t>> valB = Polynomial::valMapping(K, H, nonZeroRowsB, nonZeroColsB, p);
  Polynomial::printMapping(valB, "val_B");

  vector<vector<int64_t>> nonZeroRowsC = Polynomial::getNonZeroRows(C);
  vector<vector<int64_t>> rowC = Polynomial::createMapping(K, H, nonZeroRowsC);
  // rowC[1].push_back(1);
  // rowC[1].push_back(59);
  // rowC[1].push_back(125);
  // rowC[1].push_back(1);
  // rowC[1].push_back(135);
  // rowC[1].push_back(42);
  Polynomial::printMapping(rowC, "row_C");
  vector<vector<int64_t>> nonZeroColsC = Polynomial::getNonZeroCols(C);
  vector<vector<int64_t>> colC = Polynomial::createMapping(K, H, nonZeroColsC);
  // colC[1].push_back(125);
  // colC[1].push_back(59);
  // colC[1].push_back(1);
  // colC[1].push_back(1);
  // colC[1].push_back(42);
  // colC[1].push_back(59);
  Polynomial::printMapping(colC, "col_C");
  vector<vector<int64_t>> valC = Polynomial::valMapping(K, H, nonZeroRowsC, nonZeroColsC, p);
  Polynomial::printMapping(valC, "val_C");


  vector<int64_t> rowA_x = Polynomial::setupLagrangePolynomial(rowA[0], rowA[1], p, "rowA(x)");
  vector<int64_t> colA_x = Polynomial::setupLagrangePolynomial(colA[0], colA[1], p, "colA(x)");
  vector<int64_t> valA_x = Polynomial::setupLagrangePolynomial(valA[0], valA[1], p, "valA(x)");

  vector<int64_t> rowB_x = Polynomial::setupLagrangePolynomial(rowB[0], rowB[1], p, "rowB(x)");
  vector<int64_t> colB_x = Polynomial::setupLagrangePolynomial(colB[0], colB[1], p, "colB(x)");
  vector<int64_t> valB_x = Polynomial::setupLagrangePolynomial(valB[0], valB[1], p, "valB(x)");

  vector<int64_t> rowC_x = Polynomial::setupLagrangePolynomial(rowC[0], rowC[1], p, "rowC(x)");
  vector<int64_t> colC_x = Polynomial::setupLagrangePolynomial(colC[0], colC[1], p, "colC(x)");
  vector<int64_t> valC_x = Polynomial::setupLagrangePolynomial(valC[0], valC[1], p, "valC(x)");

  vector<int64_t> O_AHP;

  O_AHP.insert(O_AHP.end(), rowA_x.begin(), rowA_x.end());
  O_AHP.insert(O_AHP.end(), colA_x.begin(), colA_x.end());
  O_AHP.insert(O_AHP.end(), valA_x.begin(), valA_x.end());

  O_AHP.insert(O_AHP.end(), rowB_x.begin(), rowB_x.end());
  O_AHP.insert(O_AHP.end(), colB_x.begin(), colB_x.end());
  O_AHP.insert(O_AHP.end(), valB_x.begin(), valB_x.end());

  O_AHP.insert(O_AHP.end(), rowC_x.begin(), rowC_x.end());
  O_AHP.insert(O_AHP.end(), colC_x.begin(), colC_x.end());
  O_AHP.insert(O_AHP.end(), valC_x.begin(), valC_x.end());

  cout << "O_AHP = {";
  for (int64_t i = 0; i < O_AHP.size(); i++) {
    cout << O_AHP[i];
    if (i != O_AHP.size() - 1) {
      cout << ", ";
    }
  }
  cout << "}" << endl;

  int64_t Com0_AHP = 0, Com1_AHP = 0, Com2_AHP = 0, Com3_AHP = 0, Com4_AHP = 0, Com5_AHP = 0, Com6_AHP = 0, Com7_AHP = 0, Com8_AHP = 0;
  // Ensure ck and *_x vectors are of the same size before iterating
  // size_t minSize = std::min(rowA_x.size(), ck.size()); // **Changed: Added minSize for safer iteration**
  // for (int64_t i = 0; i < minSize; i++) { // **Changed: Use minSize in loop condition**
  for (int64_t i = 0; i < rowA_x.size(); i++) {
    Com0_AHP += (ck[i] * rowA_x[i]) % p;
    Com1_AHP += (ck[i] * colA_x[i]) % p;
    Com2_AHP += (ck[i] * valA_x[i]) % p;
    
    Com3_AHP += (ck[i] * rowB_x[i]) % p;
    Com4_AHP += (ck[i] * colB_x[i]) % p;
    Com5_AHP += (ck[i] * valB_x[i]) % p;
    
    Com6_AHP += (ck[i] * rowC_x[i]) % p;
    Com7_AHP += (ck[i] * colC_x[i]) % p;
    Com8_AHP += (ck[i] * valC_x[i]) % p;

    Com0_AHP %= p;
    Com1_AHP %= p;
    Com2_AHP %= p;
    Com3_AHP %= p;
    Com4_AHP %= p;
    Com5_AHP %= p;
    Com6_AHP %= p;
    Com7_AHP %= p;
    Com8_AHP %= p;
  }
  cout << "Com0_AHP = " << Com0_AHP << endl;
  cout << "Com1_AHP = " << Com1_AHP << endl;
  cout << "Com2_AHP = " << Com2_AHP << endl;
  cout << "Com3_AHP = " << Com3_AHP << endl;
  cout << "Com4_AHP = " << Com4_AHP << endl;
  cout << "Com5_AHP = " << Com5_AHP << endl;
  cout << "Com6_AHP = " << Com6_AHP << endl;
  cout << "Com7_AHP = " << Com7_AHP << endl;
  cout << "Com8_AHP = " << Com8_AHP << endl;


  // TODO: Add from device_config.json to the program_commitment.json
  /*
  device_config.json
  {
    "CommitmentID":  64-bit,
    "Class": 32-bit Integer,
    "IoT_Manufacturer_Name": String,
    "IoT_Device_Name": String,
    "Device_Hardware_Version": float,
    "Firmware_Version": float,
    "code_block": 64-bit Array
  }
// TODO: Add commitmentID to the program_commitment.json.
  
  */


  ordered_json commitment;
  commitment.clear();
  commitment["commitmentID"] = commitmentID;
  commitment["IoT_Manufacturer_Name"] = IoT_Manufacturer_Name;
  commitment["IoT_Device_Name"] = IoT_Device_Name;
  commitment["Device_Hardware_Version"] = Device_Hardware_Version;
  commitment["Firmware_Version"] = Firmware_Version;
  commitment["Class"] = Class;
  commitment["m"] = m;
  commitment["n"] = n;
  commitment["p"] = p;
  commitment["g"] = g;
  commitment["RowA"] = rowA_x;
  commitment["ColA"] = colA_x;
  commitment["ValA"] = valA_x;

  commitment["RowB"] = rowB_x;
  commitment["ColB"] = colB_x;
  commitment["ValB"] = valB_x;

  commitment["RowC"] = rowC_x;
  commitment["ColC"] = colC_x;
  commitment["ValC"] = valC_x;
  
  commitment["Curve"] = "bn128";
  commitment["PolynomialCommitment"] = "KZG";

  // Serialize JSON object to a string
  std::string commitmentString = commitment.dump();
  // Write JSON object to a file
  std::ofstream commitmentFile("data/program_commitment.json");
  if (commitmentFile.is_open()) {
      commitmentFile << commitmentString;
      commitmentFile.close();
      std::cout << "JSON data has been written to program_commitment.json\n";
  } else {
      std::cerr << "Error opening file for writing\n";
  }

  vector<vector<int64_t>> nonZeroB;
  for(int64_t i = 0; i < nonZeroRowsB[0].size(); i++){
    nonZeroB.push_back({nonZeroRowsB[0][i], nonZeroColsB[0][i], nonZeroColsB[1][i]});
  }
  ordered_json program_param;
  program_param.clear();
  program_param["Class"] = Class;
  program_param["A"] = nonZeroColsA[0];
  program_param["B"] = nonZeroB;
  program_param["rA"] = rowA[1];
  program_param["cA"] = colA[1];
  program_param["vA"] = valA[1];
  program_param["rB"] = rowB[1];
  program_param["cB"] = colB[1];
  program_param["vB"] = valB[1];
  program_param["rC"] = rowC[1];
  program_param["cC"] = colC[1];
  program_param["vC"] = valC[1];


  // Serialize JSON object to a string
  std::string program_paramString = program_param.dump();
  // Write JSON object to a file
  std::ofstream program_paramFile("data/program_param.json");
  if (program_paramFile.is_open()) {
      program_paramFile << program_paramString;
      program_paramFile.close();
      std::cout << "JSON data has been written to program_param.json\n";
  } else {
      std::cerr << "Error opening file for writing\n";
  }
}

int main() {
  std::string configFile = "device_config.json", setupFile = "data/setup3.json", assemblyFile = "program.s", newAssemblyFile = "program_new.s";

  // TODO: Remove the hard coded file names and use the inputs from user

  // std::string configFile, setupFile, assemblyFile, newAssemblyFile;
  // Input filenames
  // std::cout << "Enter the device config file name: ";
  // std::cin >> configFile;
  // std::cout << "Enter setup file name: ";
  // std::cin >> setupFile;
  // std::cout << "Enter the program assembly file name: ";
  // std::cin >> assemblyFile;
  // std::cout << "Enter the output file name for modified assembly: ";
  // std::cin >> newAssemblyFile;

  nlohmann::json config;
  auto [startLine, endLine] = parseDeviceConfig(configFile, config);

  modifyAndSaveAssembly(assemblyFile, newAssemblyFile, startLine, endLine);

  std::cout << "Modified assembly file saved as: " << newAssemblyFile << std::endl;

  // TODO: update this part to be dynamic
  commitmentGenerator(setupFile);
  return 0;
}


// Strore Matrix A and B in the program_param.json
// {
//   "A": [col1,col2,col3,...],
//   "B": [[row,col,val],[row,col,val],[row,col,val],...]
// }


/*
Read program.s and insert assembly macro SaveReg() before and after the instructions block which is specified by the code_block in the device_config.json file
Example code_block: 37-39
Also, add proofgen macro in assembly at the end of the last instuction and store it in program_new.s.

program.s   
36 JMP            
37 Add        
38 Mul        
39 Add        
40 BEQ      
41 BEG      
42 Mul      
43 SHL   

program_new.s
36 JMP 
37 saveRegX()
38 Add        
37 saveRegW0()
39 Mul        
37 saveRegW1()
40 Add
41 saveRegY()
42 proofGen()
43 BEQ      
44 BEG      
45 Mul      
46 SHL   

            
*/



