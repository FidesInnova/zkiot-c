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
mul s1, s2, s3
sw x9, x9_array(4)
addi s2, s2, 11
sw x18, x18_array(4)
mul s1, s2, s4
sw x9, x9_array(8)
mul s4, s3, s4
sw x20, x20_array(4)
jal proofGenerator

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
x0_array:    .space 4   # Array for x0
x1_array:    .space 4   # Array for x1
x2_array:    .space 4   # Array for x2
x3_array:    .space 4   # Array for x3
x4_array:    .space 4   # Array for x4
x5_array:    .space 4   # Array for x5
x6_array:    .space 4   # Array for x6
x7_array:    .space 4   # Array for x7
x8_array:    .space 4   # Array for x8
x9_array:    .space 12   # Array for x9
x10_array:    .space 4   # Array for x10
x11_array:    .space 4   # Array for x11
x12_array:    .space 4   # Array for x12
x13_array:    .space 4   # Array for x13
x14_array:    .space 4   # Array for x14
x15_array:    .space 4   # Array for x15
x16_array:    .space 4   # Array for x16
x17_array:    .space 4   # Array for x17
x18_array:    .space 8   # Array for x18
x19_array:    .space 4   # Array for x19
x20_array:    .space 8   # Array for x20
x21_array:    .space 4   # Array for x21
x22_array:    .space 4   # Array for x22
x23_array:    .space 4   # Array for x23
x24_array:    .space 4   # Array for x24
x25_array:    .space 4   # Array for x25
x26_array:    .space 4   # Array for x26
x27_array:    .space 4   # Array for x27
x28_array:    .space 4   # Array for x28
x29_array:    .space 4   # Array for x29
x30_array:    .space 4   # Array for x30
x31_array:    .space 4   # Array for x31

