#!/bin/bash

# Step 1: Compile the program.cpp to assembly
riscv64-unknown-elf-g++ -S program.cpp -o program.s -lstdc++
if [ $? -ne 0 ]; then
    echo "Compilation failed"
    exit 1
fi

# Step 2: Find the line number that uses 'mul' after '#APP' in program.s
line_number=$(awk '/#APP/{flag=1; next} flag && /mul/{print NR; exit}' program.s)
if [ -z "$line_number" ]; then
    echo "No 'mul' instruction found after '#APP'"
    exit 1
fi

# Step 3: Read the Class value from device_config.json
config_file="/root/zkiot_c_source/gcc/V2.0/device_config.json"
class_value=$(jq -r '.Class' "$config_file")
if [ -z "$class_value" ]; then
    echo "Class value not found in device_config.json"
    exit 1
fi

# Step 4: Read n_g from class.json based on the Class value
class_file="/root/zkiot_c_source/gcc/V2.0/class.json"
n_g=$(jq --arg class_value "$class_value" '.[$class_value].n_g' "$class_file")
if [ -z "$n_g" ]; then
    echo "n_g not found for Class $class_value in class.json"
    exit 1
fi

# Step 5: Calculate the new values for code_block
second_value=$((line_number + n_g -1))

# Step 6: Update device_config.json with the new values
temp_file=$(mktemp)
jq --argjson line_number "$line_number" --argjson second_value "$second_value" '.code_block = [$line_number, $second_value]' "$config_file" > "$temp_file" && mv "$temp_file" "$config_file"

# Step 7: Run the commitmentGenerator and store the output logs
# Check if the log directory exists, if not, create it
log_dir="log"
if [ ! -d "$log_dir" ]; then
    mkdir -p "$log_dir"
fi
./commitmentGenerator > log/commitmentGenerator.log 2>&1
if [ $? -ne 0 ]; then
    echo "commitmentGenerator execution failed"
    exit 1
fi

# Step 8: Build the new file using the updated code_block and store the output logs
riscv64-unknown-elf-g++ program_new.s polynomial.cpp -o program -lstdc++
if [ $? -ne 0 ]; then
    echo "Build failed"
    exit 1
fi

# Step 9: Execute the program using qemu-riscv64-static and store the output logs
qemu-riscv64-static program <<EOF > log/proofGeneration.log 2>&1
$(cat data/program_commitment.json)

$(cat data/program_param.json)

$(cat class.json)

$(cat data/setup$class_value.json)
EOF

if [ $? -ne 0 ]; then
    echo "Program execution failed"
    exit 1
fi

# Step 10: Extract the line that starts with {"commitmentID": from proofGeneration.log and store it as a JSON in data/proof.json
grep '{"commitmentID":' log/proofGeneration.log | jq . > data/proof.json

if [ $? -ne 0 ]; then
    echo "Failed to extract commitmentID JSON from proofGeneration.log"
    exit 1
fi

# Step 11: Run the verifier and store the output logs
./verifier > log/verifier.log 2>&1
if [ $? -ne 0 ]; then
    echo "Verifier execution failed"
    exit 1
fi

# Step 12: Check verifier.log for verification result
if grep -q 'verify!' log/verifier.log; then
    echo "Verification: true"
else
    echo "Verification: false"
fi

# echo "Script completed successfully"