	.file	"program.cpp"
	.option nopic
	.attribute arch, "rv32i2p1_m2p0_a2p1_f2p2_d2p2_c2p0_zicsr2p0_zifencei2p0"
	.attribute unaligned_access, 0
	.attribute stack_align, 16
	.text
	.align	1
	.globl	main
	.type	main, @function
main:
.LFB0:
	.cfi_startproc
	addi	sp,sp,-16
	.cfi_def_cfa_offset 16
	sw	s0,12(sp)
	.cfi_offset 8, -4
	addi	s0,sp,16
	.cfi_def_cfa 8, 0
 #APP
# 6 "program.cpp" 1
	li s2, 4
	li s3, 5
	li s4, 26
jal store_register_instances
mul s2, s2, 5
jal store_register_instances
addi s2, s2, 11
jal store_register_instances
mul s2, s2, 26
jal store_register_instances
addi s2, s2, 0
jal store_register_instances

# 0 "" 2
 #NO_APP
	li	a5,0
	mv	a0,a5
	lw	s0,12(sp)
	.cfi_restore 8
	.cfi_def_cfa 2, 16
	addi	sp,sp,16
	.cfi_def_cfa_offset 0
	jr	ra
	.cfi_endproc
.LFE0:
	.size	main, .-main
	.ident	"GCC: (13.2.0-11ubuntu1+12) 13.2.0"

  #### Subroutine Code (`store_registers.s`)

  ```assembly
          .data
  a0_saved:    .word 0               # Temporary storage for the original value of a0
  last_space_instance:  .word 0      # Temporary storage for the latest instance value

  x0_array:    .space 12   # Array for x0
  x1_array:    .space 12   # Array for x1
  x2_array:    .space 12   # Array for x2
  x3_array:    .space 12   # Array for x3
  x4_array:    .space 12   # Array for x4
  x5_array:    .space 12   # Array for x5
  x6_array:    .space 12   # Array for x6
  x7_array:    .space 12   # Array for x7
  x8_array:    .space 12   # Array for x8
  x9_array:    .space 12   # Array for x9
  x10_array:   .space 12   # Array for x10
  x11_array:   .space 12   # Array for x11
  x12_array:   .space 12   # Array for x12
  x13_array:   .space 12   # Array for x13
  x14_array:   .space 12   # Array for x14
  x15_array:   .space 12   # Array for x15
  x16_array:   .space 12   # Array for x16
  x17_array:   .space 12   # Array for x17
  x18_array:   .space 12   # Array for x18
  x19_array:   .space 12   # Array for x19
  x20_array:   .space 12   # Array for x20
  x21_array:   .space 12   # Array for x21
  x22_array:   .space 12   # Array for x22
  x23_array:   .space 12   # Array for x23
  x24_array:   .space 12   # Array for x24
  x25_array:   .space 12   # Array for x25
  x26_array:   .space 12   # Array for x26
  x27_array:   .space 12   # Array for x27
  x28_array:   .space 12   # Array for x28
  x29_array:   .space 12   # Array for x29
  x30_array:   .space 12   # Array for x30
  x31_array:   .space 12   # Array for x31

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
      
      sw a0, last_space_instance     # Save the original value of last instance

      # Restore original value of a0 from saved location
      lw a0, a0_saved                # Restore original value of a0

      ret                            # Return from function
  
