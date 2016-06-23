# mark_description "Intel(R) Fortran Intel(R) 64 Compiler XE for applications running on Intel(R) 64, Version 15.0.1.133 Build 2";
# mark_description "0141023";
# mark_description "-I/opt/cray/mpt/7.3.1/gni/sma/include -I/opt/cray/libsci/13.3.0/INTEL/14.0/x86_64/include -I/opt/cray/ga/5.3";
# mark_description ".0.5/INTEL/14.0/include -I/opt/cray/mpt/7.3.1/gni/mpich-intel/14.0/include -I/opt/cray/rca/1.0.0-2.0502.5721";
# mark_description "2.2.56.ari/include -I/opt/cray/alps/5.2.3-2.0502.9295.14.14.ari/include -I/opt/cray/xpmem/0.1-2.0502.57015.1";
# mark_description ".15.ari/include -I/opt/cray/gni-headers/4.0-1.0502.10317.9.2.ari/include -I/opt/cray/dmapp/7.0.1-1.0502.1024";
# mark_description "6.8.47.ari/include -I/opt/cray/pmi/5.0.10-1.0000.11050.0.0.ari/include -I/opt/cray/ugni/6.0-1.0502.10245.9.9";
# mark_description ".ari/include -I/opt/cray/udreg/2.3.2-1.0502.9889.2.20.ari/include -I/usr/local/include -I/opt/cray/wlm_detec";
# mark_description "t/1.0-1.0502.57063.1.1.ari/include -I/opt/cray/krca/1.0.0-2.0502.57202.2.45.ari/include -I/opt/cray-hss-deve";
# mark_description "l/7.2.0/include -mavx -static -D__CRAYXC -D__CRAY_SANDYBRIDGE -D__CRAYXT_COMPUTE_LINUX_TARGET --save-temps -";
# mark_description "c -o /global/project/projectdirs/m1778/malbon/surfgen/surfgen-edison/source/getver.o -DSGENVER=\"2.8.8\"";
	.file "getver.F90"
	.text
..TXTST0:
# -- Begin  getver_
	.text
# mark_begin;
       .align    16,0x90
	.globl getver_
getver_:
# parameter 1: %rdi
# parameter 2: %rsi
..B1.1:                         # Preds ..B1.0
..___tag_value_getver_.1:                                       #8.12
        movl      $72, %esi                                     #12.9
        movl      $__STRLITPACK_0, %edx                         #12.9
        movl      $5, %ecx                                      #12.9
        xorl      %r8d, %r8d                                    #12.9
        jmp       for_cpystr                                    #12.9
        .align    16,0x90
..___tag_value_getver_.3:                                       #
                                # LOE
# mark_end;
	.type	getver_,@function
	.size	getver_,.-getver_
	.data
# -- End  getver_
	.section .rodata.str1.4, "aMS",@progbits,1
	.align 4
__STRLITPACK_0:
	.long	775433778
	.word	56
	.type	__STRLITPACK_0,@object
	.size	__STRLITPACK_0,6
	.data
	.section .note.GNU-stack, ""
// -- Begin DWARF2 SEGMENT .eh_frame
	.section .eh_frame,"a",@progbits
.eh_frame_seg:
	.align 8
	.4byte 0x00000014
	.8byte 0x7801000100000000
	.8byte 0x0000019008070c10
	.4byte 0x00000000
	.4byte 0x00000014
	.4byte 0x0000001c
	.8byte ..___tag_value_getver_.1
	.8byte ..___tag_value_getver_.3-..___tag_value_getver_.1
# End
