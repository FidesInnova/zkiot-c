
// Class Fp

// function Setup()
// TargetOpcodes = (170, 234, 459, 5933, 9456)
// a0 = len (TargetOpcodes)
// l = 54
// g
// d
// pp
// Fp, example p= 181, Final Fp=2^31-2^27+1=15*2^27+1 


// Function Commit()


// Function Proof-Eval-Generation


// Function Proof-Eval-Verification


//=====Using zk-IoT-C lib=======
// Function IoT-Programmer() {
//   100,000 opocde in LLVM file => select 10 opcode
//   setup()
//   commit()
//   submit d,g,pp,commit to blockchain
//   write in a file for the IoT device
// }

// Function IoT-Device-Proof-Generator() {
//   Read d,g,pp, commit from the IoT device memory 
//   Proof-Eval-Generation <- important from computing power perspective
//   Submit result+proof to cold-db
// } => 100,000 OPcode
// (Execution 10 opcode -> 5,000,000 opcode -> STARK)
// (Execution 10 opcode -> ? opcode -> zk-IoT)

// Function Solidity-EndUSer-Proof-Verification() {
//   Read d,g,pp, commit from blockchain 
//   Read result+proof to cold-db
//   Proof-Eval-Verification
// }
