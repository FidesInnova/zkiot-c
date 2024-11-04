.data               # Data section
X: .word 0         # Reserve space for X
W0: .word 0        # Reserve space for W0
W1: .word 0        # Reserve space for W1
Y: .word 0         # Reserve space for Y
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
	li a0, 4
sw a0, X
addi a0, a0, 5
sw a0, W0
mul a0, a0, 11
sw a0, W1
addi a0, a0, 26
sw a0, Y

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
