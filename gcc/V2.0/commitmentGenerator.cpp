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


Read the Class and Lines from device_config.json to be used in the code and pass all parameters to the program_commitment.json.
  device_config.json
  {
    "Class": 32-bit Integer,
    "IoT_Manufacturer_Name": String,
    "IoT_Device_Name": String,
    "Device_Hardware_Version": float,
    "Firmware_Version": float,
    "Lines": 64-bit Array
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

using namespace std;

std::vector<std::string> instructions;
int64_t Class;

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

  int64_t startLine = config["Lines"][0].get<int64_t>();
  int64_t endLine = config["Lines"][1].get<int64_t>();
  Class = config["Class"].get<int64_t>();

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

  while (std::getline(assemblyFileStream, line)) {
    // Insert variables before the specified lines
    if (currentLineNumber == startLine) {
      newAssemblyFileStream << "jal store_register_instances\n";
      instructions.push_back(line);
    }
    else if(currentLineNumber > startLine && currentLineNumber <= endLine){
      newAssemblyFileStream << "jal store_register_instances\n";
      instructions.push_back(line);
    }
    else if (currentLineNumber == endLine + 1){
      newAssemblyFileStream << "jal store_register_instances\n";
      newAssemblyFileStream << "jal proofGenerator\n";
    }

    newAssemblyFileStream << line << std::endl;
    ++currentLineNumber;
  }

  int64_t spaceSize = (((endLine - startLine) + 2) * 4);
  std::string assemblyCode = R"(
  #### Subroutine Code (`store_registers.s`)

  ```assembly
          .data
  a0_saved:    .word 0               # Temporary storage for the original value of a0
  last_space_instance:  .word 0      # Temporary storage for the latest instance value

  x0_array:    .space {SPACE_SIZE}   # Array for x0
  x1_array:    .space {SPACE_SIZE}   # Array for x1
  x2_array:    .space {SPACE_SIZE}   # Array for x2
  x3_array:    .space {SPACE_SIZE}   # Array for x3
  x4_array:    .space {SPACE_SIZE}   # Array for x4
  x5_array:    .space {SPACE_SIZE}   # Array for x5
  x6_array:    .space {SPACE_SIZE}   # Array for x6
  x7_array:    .space {SPACE_SIZE}   # Array for x7
  x8_array:    .space {SPACE_SIZE}   # Array for x8
  x9_array:    .space {SPACE_SIZE}   # Array for x9
  x10_array:   .space {SPACE_SIZE}   # Array for x10
  x11_array:   .space {SPACE_SIZE}   # Array for x11
  x12_array:   .space {SPACE_SIZE}   # Array for x12
  x13_array:   .space {SPACE_SIZE}   # Array for x13
  x14_array:   .space {SPACE_SIZE}   # Array for x14
  x15_array:   .space {SPACE_SIZE}   # Array for x15
  x16_array:   .space {SPACE_SIZE}   # Array for x16
  x17_array:   .space {SPACE_SIZE}   # Array for x17
  x18_array:   .space {SPACE_SIZE}   # Array for x18
  x19_array:   .space {SPACE_SIZE}   # Array for x19
  x20_array:   .space {SPACE_SIZE}   # Array for x20
  x21_array:   .space {SPACE_SIZE}   # Array for x21
  x22_array:   .space {SPACE_SIZE}   # Array for x22
  x23_array:   .space {SPACE_SIZE}   # Array for x23
  x24_array:   .space {SPACE_SIZE}   # Array for x24
  x25_array:   .space {SPACE_SIZE}   # Array for x25
  x26_array:   .space {SPACE_SIZE}   # Array for x26
  x27_array:   .space {SPACE_SIZE}   # Array for x27
  x28_array:   .space {SPACE_SIZE}   # Array for x28
  x29_array:   .space {SPACE_SIZE}   # Array for x29
  x30_array:   .space {SPACE_SIZE}   # Array for x30
  x31_array:   .space {SPACE_SIZE}   # Array for x31

      .text
      .globl store_register_instances
  store_register_instances:
      sw a0, a0_saved                # Save the original value of a0

      # Load the instance index into a0 from memory
      lw a0, last_space_instance                # Read saved value of last instance

      # Store each register's value in its respective array
      sw x0, x0_array(a0)            # Store x0 in x0_array at index given by a0
      sw x1, x1_array(a0)            # Store x1 in x1_array at index given by a0
      sw x2, x2_array(a0)            # Store x2 in x2_array at index given by a0
      sw x3, x3_array(a0)            # Store x3 in x3_array at index given by a0
      sw x4, x4_array(a0)            # Store x4 in x4_array at index given by a0
      sw x5, x5_array(a0)            # Store x5 in x5_array at index given by a0
      sw x6, x6_array(a0)            # Store x6 in x6_array at index given by a0
      sw x7, x7_array(a0)            # Store x7 in x7_array at index given by a0
      sw x8, x8_array(a0)            # Store x8 in x8_array at index given by a0
      sw x9, x9_array(a0)            # Store x9 in x9_array at index given by a0
      sw x10, x10_array(a0)          # Store x10 in x10_array at index given by a0
      sw x11, x11_array(a0)          # Store x11 in x11_array at index given by a0
      sw x12, x12_array(a0)          # Store x12 in x12_array at index given by a0
      sw x13, x13_array(a0)          # Store x13 in x13_array at index given by a0
      sw x14, x14_array(a0)          # Store x14 in x14_array at index given by a0
      sw x15, x15_array(a0)          # Store x15 in x15_array at index given by a0
      sw x16, x16_array(a0)          # Store x16 in x16_array at index given by a0
      sw x17, x17_array(a0)          # Store x17 in x17_array at index given by a0
      sw x18, x18_array(a0)          # Store x18 in x18_array at index given by a0
      sw x19, x19_array(a0)          # Store x19 in x19_array at index given by a0
      sw x20, x20_array(a0)          # Store x20 in x20 array at index given by a0
      sw x21, x21_array(a0)          # Store x21 in x21_array at index given by a0
      sw x22, x22_array(a0)          # Store x22 in x22_array at index given by a0
      sw x23, x23_array(a0)          # Store x23 in x23_array at index given by a0
      sw x24, x24_array(a0)          # Store x24 in x24_array at index given by a0
      sw x25, x25_array(a0)          # Store x25 in x25_array at index given by a0
      sw x26, x26_array(a0)          # Store x26 in x26_array at index given by a0
      sw x27, x27_array(a0)          # Store x27 in x27_array at index given by a0
      sw x28, x28_array(a0)          # Store x28 in x28_array at index given by a0
      sw x29, x29_array(a0)          # Store x29 in x29_array at index given by a0
      sw x30, x30_array(a0)          # Store x30 in x30_array at index given by a0
      sw x31, x31_array(a0)          # Store x31 in x31_array at index given by a0
      
      addi a0, a0, 4
      sw a0, last_space_instance     # Save the original value of last instance

      # Restore original value of a0 from saved location
      lw a0, a0_saved                # Restore original value of a0

      ret                            # Return from function
  )";

  // Convert spaceSize to a string
  std::ostringstream oss;
  oss << spaceSize;
  std::string spaceSizeStr = oss.str();

  // Replace all instances of "{SPACE_SIZE}" with the actual value of spaceSize
  size_t pos = 0;
  while ((pos = assemblyCode.find("{SPACE_SIZE}", pos)) != std::string::npos) {
      assemblyCode.replace(pos, spaceSize, spaceSizeStr); // spaceSize is the length of "{SPACE_SIZE}"
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
    std::string operation, operationStr, leftStr, rightStr;
    
    ss >> operation >> operationStr;
    ss >> leftStr >> rightStr;
    leftStr = Polynomial::trim(leftStr);
    rightStr = Polynomial::trim(rightStr);
    leftStr = Polynomial::removeCommas(leftStr);
    rightStr = Polynomial::removeCommas(rightStr);
    // std::cout << "operation: " << operation << "\tleftInt: " << leftStr << "\trightInt: " << rightStr << "\n";
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
    std::string operation, operationStr, leftStr, rightStr;
    ss >> operation >> operationStr;
    int64_t R = std::stoi(operationStr.substr(1));

    if (operation == "add" || operation == "addi" || operation == "mul") {
      ss >> leftStr >> rightStr;

      // Remove commas
      leftStr = Polynomial::removeCommas(leftStr);
      rightStr = Polynomial::removeCommas(rightStr);
      // Trim spaces
      leftStr = Polynomial::trim(leftStr);
      rightStr = Polynomial::trim(rightStr);

      int64_t leftInt,rightInt;
      ss >> leftStr >> rightStr;
      
      C[1+n_i+i][1+n_i+i] = 1;

      if (operation == "add" || operation == "addi") {
        A[1+n_i+i][0] = 1;
        if (std::isdigit(leftStr[0])) {
          leftInt = std::stoi(leftStr);
          B[1+n_i+i][0] = leftInt;
        }
        else {
          B[1+n_i+i][1+i] = 1;
        }
        if(std::isdigit(rightStr[0])){
          rightInt = std::stoi(rightStr);
          B[1+n_i+i][0] = rightInt;
        }
        else {
          B[1+n_i+i][1+i] = 1;
        }

    } else if (operation == "mul") {
        if (std::isdigit(leftStr[0])) {
          leftInt = std::stoi(leftStr);
          A[1+n_i+i][0] = leftInt;
        }
        else {
          A[1+n_i+i][1+i] = 1;
        }
        if (std::isdigit(rightStr[0])) {
          rightInt = std::stoi(rightStr);
          B[1+n_i+i][0] = rightInt;
        }
        else {
          B[1+n_i+i][1+i] = 1;
        }
      }
    }
    
    else {
      cout << "!!! Undefined instruction in the defiend Line range !!!\n" << operation << endl;
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
    "Lines": 64-bit Array
  }
// TODO: Add commitmentID to the program_commitment.json.
  
  */


  ordered_json commitment;
  commitment.clear();
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
  std::string configFile = "device_config.json", setupFile = "setup3.json", assemblyFile = "program.s", newAssemblyFile = "new_program.s";

  // TODO: Remove the hard coded file names and use the inputs from user

  // Input filenames
  // std::cout << "Enter the device config file name: ";
  // std::cin >> configFile;
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
Read program.s and insert assembly macro SaveReg() before and after the instructions block which is specified by the Lines in the device_config.json file
Example Lines: 37-39
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



