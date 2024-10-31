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

#include "fidesinnova.h"

void fidesinnova::commitmentGenerator(String path, int64_t g, int64_t p) {
  //read /setup.json file to initialize the cryptographic environment
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

    //I add this for error handling for if there is no ''
    // if (firstSpace == -1) {
    //   Serial.println("Error: Invalid instruction format");
    //   return;  // Exit if format is wrong
    // }

    String operation = inst.substring(0, firstSpace);   // Extract the operation
    String remainder = inst.substring(firstSpace + 1);  // Extract the remainder (registers and values)
    
    // //add .trim() to achieve a better format of input
    // String remainder = inst.substring(firstSpace + 1).trim();  // Extract the remainder (registers and values)

    int64_t secondSpace = remainder.indexOf(' ');
    String rStr = remainder.substring(0, secondSpace);   // Extract first register or value
    String rest = remainder.substring(secondSpace + 1);  // The rest of the instruction

    //String rest = remainder.substring(secondSpace + 1).trim();

    if (operation == "li") {
      n_i++;
      String xStr = rest;  // Only one value left for li

      int64_t R = rStr.substring(1).toInt();  // Get the register number

      //I add this to chekch if it has been successfully converted to int, if not, R = 0;
      // if (R == 0 && rStr.substring(1) != "0") {
      //   Serial.println("Error: Invalid register number");
      //   return;
      // }

      int64_t X = xStr.toInt();               // Get the immediate value

      // if (X == 0 && xStr != "0") {
      //   Serial.println("Error: Invalid immediate value");
      //   return;
      // }

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
  m = (((Polynomial::power(n, 2, p) - n) / 2) - ((Polynomial::power(t, 2, p) - t) / 2)) % p;

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

    //constrct A,B,C Matrix
    if (operation == "li") {
      X = xStr.toInt();
      z[i + 1] = X % p;
    } else {
      gateIndex++;
      int64_t newI = gateIndex + n_i;
      C[newI][newI] = 1;

      if (operation == "addi") {
        A[n_i + newI - 1][0] = 1;

        if (isDigit(xStr[0])) {
          X = xStr.toInt();
          z[i + 1] = (z[i] + X) % p;
          B[n_i + newI - 1][0] = X;
          B[n_i + newI - 1][newI - 1] = 1;
        } else if (isDigit(yStr[0])) {
          Y = yStr.toInt();
          z[i + 1] = (z[i] + Y) % p;
          B[n_i + newI - 1][0] = Y;
          B[n_i + newI - 1][newI - 1] = 1;
        }
      } else if (operation == "mul") {
        if (isDigit(xStr[0])) {
          X = xStr.toInt();
          z[i + 1] = (z[i] * X) % p;
          A[n_i + newI - 1][newI - 1] = X;
          B[n_i + newI - 1][newI - 1] = 1;
        } else if (isDigit(yStr[0])) {
          Y = yStr.toInt();
          z[i + 1] = (z[i] * Y) % p;
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

  // Vector H to store powers of w
  vector<int64_t> H;
  int64_t w, g_n;

  H.push_back(1);
  g_n = ((p - 1) / n) % p;
  w = Polynomial::power(g, g_n, p);
  for (int64_t i = 1; i < n; i++) {
    H.push_back(Polynomial::power(w, i, p));
  }
  Serial.print("H[n]: ");
  for (int64_t i = 0; i < n; i++) {
    Serial.print(H[i]);
    Serial.print(" ");
  }
  Serial.println("");

  int64_t y, g_m;

  // Vector K to store powers of y
  vector<int64_t> K;
  K.push_back(1);
  g_m = ((p - 1) * Polynomial::pInverse(m, p)) % p;
  y = Polynomial::power(g, g_m, p);
  for (int64_t i = 1; i < m; i++) {
    K.push_back(Polynomial::power(y, i, p));
  }
  Serial.print("K[m]: ");
  for (int64_t i = 0; i < m; i++) {
    Serial.print(K[i]);
    Serial.print(" ");
  }
  Serial.println("");
  
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
  vector<vector<int64_t>> valA = Polynomial::valMapping(K, H, nonZeroRowsA, nonZeroColsA, p);
  Polynomial::printMapping(valA, "val_A");

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
  vector<vector<int64_t>> valB = Polynomial::valMapping(K, H, nonZeroRowsB, nonZeroColsB, p);
  Polynomial::printMapping(valB, "val_B");

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

  Serial.print("O_AHP = {");
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

  for (int64_t i = 0; i < O_AHP.size(); i++) {
    Serial.print(O_AHP[i]);
    if (i != O_AHP.size() - 1) {
      Serial.print(", ");
    }
  }
  Serial.println("}");

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
  Serial.print("Com0_AHP = ");
  Serial.println(Com0_AHP);

  Serial.print("Com1_AHP = ");
  Serial.println(Com1_AHP);

  Serial.print("Com2_AHP = ");
  Serial.println(Com2_AHP);

  Serial.print("Com3_AHP = ");
  Serial.println(Com3_AHP);

  Serial.print("Com4_AHP = ");
  Serial.println(Com4_AHP);

  Serial.print("Com5_AHP = ");
  Serial.println(Com5_AHP);

  Serial.print("Com6_AHP = ");
  Serial.println(Com6_AHP);

  Serial.print("Com7_AHP = ");
  Serial.println(Com7_AHP);

  Serial.print("Com8_AHP = ");
  Serial.println(Com8_AHP);


  // TODO: Add from device_config.json to the program_commitment.json
  /*
  device_config.jso
  {
    "CommitmentID":  64-bit,
    "Class": 32-bit Integer,
    "IoT_Manufacturer_Name": String,
    "IoT_Device_Name": String,
    "Device_Hardware_Version": float,
    "Firmware_Version": float,
    "Lines": 64-bit Array
  }
  
  */




  DynamicJsonDocument doc(2048);
  JsonArray jsonArray;
  jsonArray = doc.createNestedArray("m");
  jsonArray.add(m);

  jsonArray = doc.createNestedArray("n");
  jsonArray.add(n);

  jsonArray = doc.createNestedArray("RowA");
  for (int64_t value : rowA_x) {
    jsonArray.add(value);
  }
  jsonArray = doc.createNestedArray("ColA");
  for (int64_t value : colA_x) {
    jsonArray.add(value);
  }
  jsonArray = doc.createNestedArray("ValA");
  for (int64_t value : valA_x) {
    jsonArray.add(value);
  }

  jsonArray = doc.createNestedArray("RowB");
  for (int64_t value : rowB_x) {
    jsonArray.add(value);
  }
  jsonArray = doc.createNestedArray("ColB");
  for (int64_t value : colB_x) {
    jsonArray.add(value);
  }
  jsonArray = doc.createNestedArray("ValB");
  for (int64_t value : valB_x) {
    jsonArray.add(value);
  }

  jsonArray = doc.createNestedArray("RowC");
  for (int64_t value : rowC_x) {
    jsonArray.add(value);
  }
  jsonArray = doc.createNestedArray("ColC");
  for (int64_t value : colC_x) {
    jsonArray.add(value);
  }
  jsonArray = doc.createNestedArray("ValC");
  for (int64_t value : valC_x) {
    jsonArray.add(value);
  }

  jsonArray = doc.createNestedArray("Curve");
  jsonArray.add("bn128");

  jsonArray = doc.createNestedArray("PolynomialCommitment");
  jsonArray.add("KZG");

  String output;
  serializeJson(doc, output);
  removeFile("/program_commitment.json");
  writeFile("/program_commitment.json", output);
}
