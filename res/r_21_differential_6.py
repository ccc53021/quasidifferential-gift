# !/user/bin/env python3
# -*- coding: utf-8 -*-

# differential 6 in this paper
# input difference is 0x00000000000000000000000000000a00
# output difference is 0x00000000000000000002200000011000
# p = 2^{-128.64} : 131 * 2 + 132 * 4 + 133 * 2 + 134 * 2 + 135 * 6 = 128.64
differential_6_21r = [
    # Solution: of weight 130
    [
        # No route with weight equal to 130.
    ],
    # Solution: of weight 131
    [
        # Solution: [#1 of weight 131 and active sbox number 4]
        [
            0x00000000000000000000000000000a00,
            0x00000000000000000000000000000100,
            0x00000000000000010000000000000000,
            0x00000000000000080000000000000000,
            0x00080000000000000000000000000000,
            0x00030000000000000000000000000000,
            0x00000000000000002000000010000000,
            0x00000000000000006000000060000000,
            0x00004040000020200000000000000000,
            0x00005050000050500000000000000000,
            0x05050000000000000505000000000000,
            0x08020000000000000802000000000000,
            0x0000000000000000a000a00000000000,
            0x00000000000000001000100000000000,
            0x00000000000000000000110000000000,
            0x00000000000000000000880000000000,
            0x00000000000000000000080000000800,
            0x00000000000000000000030000000300,
            0x00000202000001010000000000000000,
            0x00000505000005050000000000000000,
            0x00000000050500000000000005050000,
            0x00000000080200000000000008020000,
            0x000000000000000000a000a000000000,
            0x00000000000000000010001000000000,
            0x00001100000000000000000000000000,
            0x00009800000000000000000000000000,
            0x00000000000000000900000008000000,
            0x00000000000000000100000003000000,
            0x00000020000010100000000000000000,
            0x00000060000080600000000000000000,
            0x000000000000000004040000020a0000,
            0x00000000000000000505000005010000,
            0x00000000000050100000000000005050,
            0x00000000000020800000000000002080,
            0x00000000000a000a0000000000000000,
            0x00000000000100010000000000000000,
            0x00000000000000000000000000110000,
            0x00000000000000000000000000850000,
            0x00000000000000c00000000000000010,
            0x00000000000000400000000000000080,
            0x00000000000000080004000000000000,
            0x00000000000000030003000000000000,
            0x00000000000000000002200000011000
        ],
        # Solution: [#2 of weight 131 and active sbox number 4]
        [
            0x00000000000000000000000000000a00,
            0x00000000000000000000000000000100,
            0x00000000000000010000000000000000,
            0x00000000000000080000000000000000,
            0x00080000000000000000000000000000,
            0x00030000000000000000000000000000,
            0x00000000000000002000000010000000,
            0x00000000000000006000000060000000,
            0x00004040000020200000000000000000,
            0x00005050000050500000000000000000,
            0x05050000000000000505000000000000,
            0x08020000000000000802000000000000,
            0x0000000000000000a000a00000000000,
            0x00000000000000001000100000000000,
            0x00000000000000000000110000000000,
            0x00000000000000000000880000000000,
            0x00000000000000000000080000000800,
            0x00000000000000000000030000000300,
            0x00000202000001010000000000000000,
            0x00000505000005050000000000000000,
            0x00000000050500000000000005050000,
            0x00000000080200000000000008020000,
            0x000000000000000000a000a000000000,
            0x00000000000000000010001000000000,
            0x00001100000000000000000000000000,
            0x00009800000000000000000000000000,
            0x00000000000000000900000008000000,
            0x00000000000000000100000003000000,
            0x00000020000010100000000000000000,
            0x00000060000080600000000000000000,
            0x000000000000000004040000020a0000,
            0x00000000000000000505000005010000,
            0x00000000000050100000000000005050,
            0x00000000000020800000000000002080,
            0x00000000000a000a0000000000000000,
            0x00000000000100010000000000000000,
            0x00000000000000000000000000110000,
            0x00000000000000000000000000a50000,
            0x00000000000000c00000000000000030,
            0x00000000000000400000000000000080,
            0x00000000000000080004000000000000,
            0x00000000000000030003000000000000,
            0x00000000000000000002200000011000
        ],
        # No route with weight equal to 131.
    ],
    # Solution: of weight 132
    [
        # Solution: [#1 of weight 132 and active sbox number 4]
        [
            0x00000000000000000000000000000a00,
            0x00000000000000000000000000000100,
            0x00000000000000010000000000000000,
            0x00000000000000080000000000000000,
            0x00080000000000000000000000000000,
            0x00030000000000000000000000000000,
            0x00000000000000002000000010000000,
            0x00000000000000006000000060000000,
            0x00004040000020200000000000000000,
            0x00005050000050500000000000000000,
            0x05050000000000000505000000000000,
            0x08020000000000000802000000000000,
            0x0000000000000000a000a00000000000,
            0x00000000000000001000100000000000,
            0x00000000000000000000110000000000,
            0x00000000000000000000880000000000,
            0x00000000000000000000080000000800,
            0x00000000000000000000030000000300,
            0x00000202000001010000000000000000,
            0x00000505000005050000000000000000,
            0x00000000050500000000000005050000,
            0x00000000080200000000000008020000,
            0x000000000000000000a000a000000000,
            0x00000000000000000010001000000000,
            0x00001100000000000000000000000000,
            0x00009800000000000000000000000000,
            0x00000000000000000900000008000000,
            0x00000000000000000100000003000000,
            0x00000020000010100000000000000000,
            0x00000060000090600000000000000000,
            0x000000000000000004050000020a0000,
            0x00000000000000000505000005010000,
            0x00000000000050100000000000005050,
            0x00000000000020800000000000002080,
            0x00000000000a000a0000000000000000,
            0x00000000000100010000000000000000,
            0x00000000000000000000000000110000,
            0x00000000000000000000000000a50000,
            0x00000000000000c00000000000000030,
            0x00000000000000400000000000000080,
            0x00000000000000080004000000000000,
            0x00000000000000030003000000000000,
            0x00000000000000000002200000011000
        ],
        # Solution: [#2 of weight 132 and active sbox number 4]
        [
            0x00000000000000000000000000000a00,
            0x00000000000000000000000000000100,
            0x00000000000000010000000000000000,
            0x00000000000000080000000000000000,
            0x00080000000000000000000000000000,
            0x00030000000000000000000000000000,
            0x00000000000000002000000010000000,
            0x00000000000000006000000060000000,
            0x00004040000020200000000000000000,
            0x00005050000050500000000000000000,
            0x05050000000000000505000000000000,
            0x08020000000000000802000000000000,
            0x0000000000000000a000a00000000000,
            0x00000000000000001000100000000000,
            0x00000000000000000000110000000000,
            0x00000000000000000000980000000000,
            0x00000000000000000000090000000800,
            0x00000000000000000000030000000300,
            0x00000202000001010000000000000000,
            0x00000505000005050000000000000000,
            0x00000000050500000000000005050000,
            0x00000000080200000000000008020000,
            0x000000000000000000a000a000000000,
            0x00000000000000000010001000000000,
            0x00001100000000000000000000000000,
            0x00009800000000000000000000000000,
            0x00000000000000000900000008000000,
            0x00000000000000000100000003000000,
            0x00000020000010100000000000000000,
            0x00000060000080600000000000000000,
            0x000000000000000004040000020a0000,
            0x00000000000000000505000005010000,
            0x00000000000050100000000000005050,
            0x00000000000020800000000000002080,
            0x00000000000a000a0000000000000000,
            0x00000000000100010000000000000000,
            0x00000000000000000000000000110000,
            0x00000000000000000000000000a50000,
            0x00000000000000c00000000000000030,
            0x00000000000000400000000000000080,
            0x00000000000000080004000000000000,
            0x00000000000000030003000000000000,
            0x00000000000000000002200000011000
        ],
        # Solution: [#3 of weight 132 and active sbox number 4]
        [
            0x00000000000000000000000000000a00,
            0x00000000000000000000000000000100,
            0x00000000000000010000000000000000,
            0x00000000000000080000000000000000,
            0x00080000000000000000000000000000,
            0x00030000000000000000000000000000,
            0x00000000000000002000000010000000,
            0x00000000000000006000000060000000,
            0x00004040000020200000000000000000,
            0x00005050000050500000000000000000,
            0x05050000000000000505000000000000,
            0x08020000000000000802000000000000,
            0x0000000000000000a000a00000000000,
            0x00000000000000001000100000000000,
            0x00000000000000000000110000000000,
            0x00000000000000000000980000000000,
            0x00000000000000000000090000000800,
            0x00000000000000000000030000000300,
            0x00000202000001010000000000000000,
            0x00000505000005050000000000000000,
            0x00000000050500000000000005050000,
            0x00000000080200000000000008020000,
            0x000000000000000000a000a000000000,
            0x00000000000000000010001000000000,
            0x00001100000000000000000000000000,
            0x00009800000000000000000000000000,
            0x00000000000000000900000008000000,
            0x00000000000000000100000003000000,
            0x00000020000010100000000000000000,
            0x00000060000080600000000000000000,
            0x000000000000000004040000020a0000,
            0x00000000000000000505000005010000,
            0x00000000000050100000000000005050,
            0x00000000000020800000000000002080,
            0x00000000000a000a0000000000000000,
            0x00000000000100010000000000000000,
            0x00000000000000000000000000110000,
            0x00000000000000000000000000850000,
            0x00000000000000c00000000000000010,
            0x00000000000000400000000000000080,
            0x00000000000000080004000000000000,
            0x00000000000000030003000000000000,
            0x00000000000000000002200000011000
        ],
        # Solution: [#4 of weight 132 and active sbox number 4]
        [
            0x00000000000000000000000000000a00,
            0x00000000000000000000000000000100,
            0x00000000000000010000000000000000,
            0x00000000000000080000000000000000,
            0x00080000000000000000000000000000,
            0x00030000000000000000000000000000,
            0x00000000000000002000000010000000,
            0x00000000000000006000000060000000,
            0x00004040000020200000000000000000,
            0x00005050000050500000000000000000,
            0x05050000000000000505000000000000,
            0x08020000000000000802000000000000,
            0x0000000000000000a000a00000000000,
            0x00000000000000001000100000000000,
            0x00000000000000000000110000000000,
            0x00000000000000000000880000000000,
            0x00000000000000000000080000000800,
            0x00000000000000000000030000000300,
            0x00000202000001010000000000000000,
            0x00000505000005050000000000000000,
            0x00000000050500000000000005050000,
            0x00000000080200000000000008020000,
            0x000000000000000000a000a000000000,
            0x00000000000000000010001000000000,
            0x00001100000000000000000000000000,
            0x00009800000000000000000000000000,
            0x00000000000000000900000008000000,
            0x00000000000000000100000003000000,
            0x00000020000010100000000000000000,
            0x00000060000090600000000000000000,
            0x000000000000000004050000020a0000,
            0x00000000000000000505000005010000,
            0x00000000000050100000000000005050,
            0x00000000000020800000000000002080,
            0x00000000000a000a0000000000000000,
            0x00000000000100010000000000000000,
            0x00000000000000000000000000110000,
            0x00000000000000000000000000850000,
            0x00000000000000c00000000000000010,
            0x00000000000000400000000000000080,
            0x00000000000000080004000000000000,
            0x00000000000000030003000000000000,
            0x00000000000000000002200000011000
        ],
        # No route with weight equal to 132.
    ],
    # Solution: of weight 133
    [
        # Solution: [#1 of weight 133 and active sbox number 4]
        [
            0x00000000000000000000000000000a00,
            0x00000000000000000000000000000100,
            0x00000000000000010000000000000000,
            0x00000000000000080000000000000000,
            0x00080000000000000000000000000000,
            0x00030000000000000000000000000000,
            0x00000000000000002000000010000000,
            0x00000000000000006000000060000000,
            0x00004040000020200000000000000000,
            0x00005050000050500000000000000000,
            0x05050000000000000505000000000000,
            0x08020000000000000802000000000000,
            0x0000000000000000a000a00000000000,
            0x00000000000000001000100000000000,
            0x00000000000000000000110000000000,
            0x00000000000000000000980000000000,
            0x00000000000000000000090000000800,
            0x00000000000000000000030000000300,
            0x00000202000001010000000000000000,
            0x00000505000005050000000000000000,
            0x00000000050500000000000005050000,
            0x00000000080200000000000008020000,
            0x000000000000000000a000a000000000,
            0x00000000000000000010001000000000,
            0x00001100000000000000000000000000,
            0x00009800000000000000000000000000,
            0x00000000000000000900000008000000,
            0x00000000000000000100000003000000,
            0x00000020000010100000000000000000,
            0x00000060000090600000000000000000,
            0x000000000000000004050000020a0000,
            0x00000000000000000505000005010000,
            0x00000000000050100000000000005050,
            0x00000000000020800000000000002080,
            0x00000000000a000a0000000000000000,
            0x00000000000100010000000000000000,
            0x00000000000000000000000000110000,
            0x00000000000000000000000000a50000,
            0x00000000000000c00000000000000030,
            0x00000000000000400000000000000080,
            0x00000000000000080004000000000000,
            0x00000000000000030003000000000000,
            0x00000000000000000002200000011000
        ],
        # Solution: [#2 of weight 133 and active sbox number 4]
        [
            0x00000000000000000000000000000a00,
            0x00000000000000000000000000000100,
            0x00000000000000010000000000000000,
            0x00000000000000080000000000000000,
            0x00080000000000000000000000000000,
            0x00030000000000000000000000000000,
            0x00000000000000002000000010000000,
            0x00000000000000006000000060000000,
            0x00004040000020200000000000000000,
            0x00005050000050500000000000000000,
            0x05050000000000000505000000000000,
            0x08020000000000000802000000000000,
            0x0000000000000000a000a00000000000,
            0x00000000000000001000100000000000,
            0x00000000000000000000110000000000,
            0x00000000000000000000980000000000,
            0x00000000000000000000090000000800,
            0x00000000000000000000030000000300,
            0x00000202000001010000000000000000,
            0x00000505000005050000000000000000,
            0x00000000050500000000000005050000,
            0x00000000080200000000000008020000,
            0x000000000000000000a000a000000000,
            0x00000000000000000010001000000000,
            0x00001100000000000000000000000000,
            0x00009800000000000000000000000000,
            0x00000000000000000900000008000000,
            0x00000000000000000100000003000000,
            0x00000020000010100000000000000000,
            0x00000060000090600000000000000000,
            0x000000000000000004050000020a0000,
            0x00000000000000000505000005010000,
            0x00000000000050100000000000005050,
            0x00000000000020800000000000002080,
            0x00000000000a000a0000000000000000,
            0x00000000000100010000000000000000,
            0x00000000000000000000000000110000,
            0x00000000000000000000000000850000,
            0x00000000000000c00000000000000010,
            0x00000000000000400000000000000080,
            0x00000000000000080004000000000000,
            0x00000000000000030003000000000000,
            0x00000000000000000002200000011000
        ],
        # No route with weight equal to 133.
    ],
    # Solution: of weight 134
    [
        # Solution: [#1 of weight 134 and active sbox number 4]
        [
            0x00000000000000000000000000000a00,
            0x00000000000000000000000000000100,
            0x00000000000000010000000000000000,
            0x00000000000000080000000000000000,
            0x00080000000000000000000000000000,
            0x00030000000000000000000000000000,
            0x00000000000000002000000010000000,
            0x00000000000000006000000060000000,
            0x00004040000020200000000000000000,
            0x00005050000050500000000000000000,
            0x05050000000000000505000000000000,
            0x08020000000000000802000000000000,
            0x0000000000000000a000a00000000000,
            0x00000000000000001000100000000000,
            0x00000000000000000000110000000000,
            0x00000000000000000000880000000000,
            0x00000000000000000000080000000800,
            0x00000000000000000000030000000300,
            0x00000202000001010000000000000000,
            0x00000505000005050000000000000000,
            0x00000000050500000000000005050000,
            0x00000000080200000000000008020000,
            0x000000000000000000a000a000000000,
            0x00000000000000000010001000000000,
            0x00001100000000000000000000000000,
            0x00008800000000000000000000000000,
            0x00000000000000000800000008000000,
            0x00000000000000000300000003000000,
            0x00002020000010100000000000000000,
            0x00009060000080600000000000000000,
            0x0000000000000000050400000a0a0000,
            0x00000000000000000505000001010000,
            0x00000000000050100000000000005010,
            0x00000000000020800000000000002080,
            0x00000000000a000a0000000000000000,
            0x00000000000100010000000000000000,
            0x00000000000000000000000000110000,
            0x00000000000000000000000000850000,
            0x00000000000000c00000000000000010,
            0x00000000000000400000000000000080,
            0x00000000000000080004000000000000,
            0x00000000000000030003000000000000,
            0x00000000000000000002200000011000
        ],
        # Solution: [#2 of weight 134 and active sbox number 4]
        [
            0x00000000000000000000000000000a00,
            0x00000000000000000000000000000100,
            0x00000000000000010000000000000000,
            0x00000000000000080000000000000000,
            0x00080000000000000000000000000000,
            0x00030000000000000000000000000000,
            0x00000000000000002000000010000000,
            0x00000000000000006000000060000000,
            0x00004040000020200000000000000000,
            0x00005050000050500000000000000000,
            0x05050000000000000505000000000000,
            0x08020000000000000802000000000000,
            0x0000000000000000a000a00000000000,
            0x00000000000000001000100000000000,
            0x00000000000000000000110000000000,
            0x00000000000000000000880000000000,
            0x00000000000000000000080000000800,
            0x00000000000000000000030000000300,
            0x00000202000001010000000000000000,
            0x00000505000005050000000000000000,
            0x00000000050500000000000005050000,
            0x00000000080200000000000008020000,
            0x000000000000000000a000a000000000,
            0x00000000000000000010001000000000,
            0x00001100000000000000000000000000,
            0x00008800000000000000000000000000,
            0x00000000000000000800000008000000,
            0x00000000000000000300000003000000,
            0x00002020000010100000000000000000,
            0x00009060000080600000000000000000,
            0x0000000000000000050400000a0a0000,
            0x00000000000000000505000001010000,
            0x00000000000050100000000000005010,
            0x00000000000020800000000000002080,
            0x00000000000a000a0000000000000000,
            0x00000000000100010000000000000000,
            0x00000000000000000000000000110000,
            0x00000000000000000000000000a50000,
            0x00000000000000c00000000000000030,
            0x00000000000000400000000000000080,
            0x00000000000000080004000000000000,
            0x00000000000000030003000000000000,
            0x00000000000000000002200000011000
        ],
        # No route with weight equal to 134.
    ],
    # Solution: of weight 135
    [
        # Solution: [#1 of weight 135 and active sbox number 4]
        [
            0x00000000000000000000000000000a00,
            0x00000000000000000000000000000100,
            0x00000000000000010000000000000000,
            0x00000000000000080000000000000000,
            0x00080000000000000000000000000000,
            0x00030000000000000000000000000000,
            0x00000000000000002000000010000000,
            0x00000000000000006000000060000000,
            0x00004040000020200000000000000000,
            0x00005050000050500000000000000000,
            0x05050000000000000505000000000000,
            0x08020000000000000802000000000000,
            0x0000000000000000a000a00000000000,
            0x00000000000000001000100000000000,
            0x00000000000000000000110000000000,
            0x00000000000000000000880000000000,
            0x00000000000000000000080000000800,
            0x00000000000000000000030000000300,
            0x00000202000001010000000000000000,
            0x00000505000005050000000000000000,
            0x00000000050500000000000005050000,
            0x00000000080200000000000008020000,
            0x000000000000000000a000a000000000,
            0x00000000000000000010001000000000,
            0x00001100000000000000000000000000,
            0x00008800000000000000000000000000,
            0x00000000000000000800000008000000,
            0x00000000000000000300000003000000,
            0x00002020000010100000000000000000,
            0x00009060000090600000000000000000,
            0x0000000000000000050500000a0a0000,
            0x00000000000000000505000001010000,
            0x00000000000050100000000000005010,
            0x00000000000020800000000000002080,
            0x00000000000a000a0000000000000000,
            0x00000000000100010000000000000000,
            0x00000000000000000000000000110000,
            0x00000000000000000000000000a50000,
            0x00000000000000c00000000000000030,
            0x00000000000000400000000000000080,
            0x00000000000000080004000000000000,
            0x00000000000000030003000000000000,
            0x00000000000000000002200000011000
        ],
        # Solution: [#2 of weight 135 and active sbox number 4]
        [
            0x00000000000000000000000000000a00,
            0x00000000000000000000000000000100,
            0x00000000000000010000000000000000,
            0x00000000000000080000000000000000,
            0x00080000000000000000000000000000,
            0x00030000000000000000000000000000,
            0x00000000000000002000000010000000,
            0x00000000000000006000000060000000,
            0x00004040000020200000000000000000,
            0x00005050000050500000000000000000,
            0x05050000000000000505000000000000,
            0x08020000000000000802000000000000,
            0x0000000000000000a000a00000000000,
            0x00000000000000001000100000000000,
            0x00000000000000000000110000000000,
            0x00000000000000000000980000000000,
            0x00000000000000000000090000000800,
            0x00000000000000000000030000000300,
            0x00000202000001010000000000000000,
            0x00000505000005050000000000000000,
            0x00000000050500000000000005050000,
            0x00000000080200000000000008020000,
            0x000000000000000000a000a000000000,
            0x00000000000000000010001000000000,
            0x00001100000000000000000000000000,
            0x00008800000000000000000000000000,
            0x00000000000000000800000008000000,
            0x00000000000000000300000003000000,
            0x00002020000010100000000000000000,
            0x00009060000080600000000000000000,
            0x0000000000000000050400000a0a0000,
            0x00000000000000000505000001010000,
            0x00000000000050100000000000005010,
            0x00000000000020800000000000002080,
            0x00000000000a000a0000000000000000,
            0x00000000000100010000000000000000,
            0x00000000000000000000000000110000,
            0x00000000000000000000000000a50000,
            0x00000000000000c00000000000000030,
            0x00000000000000400000000000000080,
            0x00000000000000080004000000000000,
            0x00000000000000030003000000000000,
            0x00000000000000000002200000011000
        ],
        # Solution: [#3 of weight 135 and active sbox number 4]
        [
            0x00000000000000000000000000000a00,
            0x00000000000000000000000000000100,
            0x00000000000000010000000000000000,
            0x00000000000000080000000000000000,
            0x00080000000000000000000000000000,
            0x00030000000000000000000000000000,
            0x00000000000000002000000010000000,
            0x00000000000000006000000060000000,
            0x00004040000020200000000000000000,
            0x00005050000050500000000000000000,
            0x05050000000000000505000000000000,
            0x08020000000000000802000000000000,
            0x0000000000000000a000a00000000000,
            0x00000000000000001000100000000000,
            0x00000000000000000000110000000000,
            0x00000000000000000000980000000000,
            0x00000000000000000000090000000800,
            0x00000000000000000000030000000300,
            0x00000202000001010000000000000000,
            0x00000505000005050000000000000000,
            0x00000000050500000000000005050000,
            0x00000000080200000000000008020000,
            0x000000000000000000a000a000000000,
            0x00000000000000000010001000000000,
            0x00001100000000000000000000000000,
            0x00008800000000000000000000000000,
            0x00000000000000000800000008000000,
            0x00000000000000000300000003000000,
            0x00002020000010100000000000000000,
            0x00009060000080600000000000000000,
            0x0000000000000000050400000a0a0000,
            0x00000000000000000505000001010000,
            0x00000000000050100000000000005010,
            0x00000000000020800000000000002080,
            0x00000000000a000a0000000000000000,
            0x00000000000100010000000000000000,
            0x00000000000000000000000000110000,
            0x00000000000000000000000000850000,
            0x00000000000000c00000000000000010,
            0x00000000000000400000000000000080,
            0x00000000000000080004000000000000,
            0x00000000000000030003000000000000,
            0x00000000000000000002200000011000
        ],
        # Solution: [#4 of weight 135 and active sbox number 4]
        [
            0x00000000000000000000000000000a00,
            0x00000000000000000000000000000100,
            0x00000000000000010000000000000000,
            0x00000000000000080000000000000000,
            0x00080000000000000000000000000000,
            0x00030000000000000000000000000000,
            0x00000000000000002000000010000000,
            0x00000000000000006000000060000000,
            0x00004040000020200000000000000000,
            0x00005050000050500000000000000000,
            0x05050000000000000505000000000000,
            0x08020000000000000802000000000000,
            0x0000000000000000a000a00000000000,
            0x00000000000000001000100000000000,
            0x00000000000000000000110000000000,
            0x00000000000000000000880000000000,
            0x00000000000000000000080000000800,
            0x00000000000000000000030000000300,
            0x00000202000001010000000000000000,
            0x00000505000005050000000000000000,
            0x00000000050500000000000005050000,
            0x00000000080200000000000008020000,
            0x000000000000000000a000a000000000,
            0x00000000000000000010001000000000,
            0x00001100000000000000000000000000,
            0x00009800000000000000000000000000,
            0x00000000000000000900000008000000,
            0x00000000000000000300000003000000,
            0x00002020000010100000000000000000,
            0x00009060000080600000000000000000,
            0x0000000000000000050400000a0a0000,
            0x00000000000000000505000001010000,
            0x00000000000050100000000000005010,
            0x00000000000020800000000000002080,
            0x00000000000a000a0000000000000000,
            0x00000000000100010000000000000000,
            0x00000000000000000000000000110000,
            0x00000000000000000000000000a50000,
            0x00000000000000c00000000000000030,
            0x00000000000000400000000000000080,
            0x00000000000000080004000000000000,
            0x00000000000000030003000000000000,
            0x00000000000000000002200000011000
        ],
        # Solution: [#5 of weight 135 and active sbox number 4]
        [
            0x00000000000000000000000000000a00,
            0x00000000000000000000000000000100,
            0x00000000000000010000000000000000,
            0x00000000000000080000000000000000,
            0x00080000000000000000000000000000,
            0x00030000000000000000000000000000,
            0x00000000000000002000000010000000,
            0x00000000000000006000000060000000,
            0x00004040000020200000000000000000,
            0x00005050000050500000000000000000,
            0x05050000000000000505000000000000,
            0x08020000000000000802000000000000,
            0x0000000000000000a000a00000000000,
            0x00000000000000001000100000000000,
            0x00000000000000000000110000000000,
            0x00000000000000000000880000000000,
            0x00000000000000000000080000000800,
            0x00000000000000000000030000000300,
            0x00000202000001010000000000000000,
            0x00000505000005050000000000000000,
            0x00000000050500000000000005050000,
            0x00000000080200000000000008020000,
            0x000000000000000000a000a000000000,
            0x00000000000000000010001000000000,
            0x00001100000000000000000000000000,
            0x00008800000000000000000000000000,
            0x00000000000000000800000008000000,
            0x00000000000000000300000003000000,
            0x00002020000010100000000000000000,
            0x00009060000090600000000000000000,
            0x0000000000000000050500000a0a0000,
            0x00000000000000000505000001010000,
            0x00000000000050100000000000005010,
            0x00000000000020800000000000002080,
            0x00000000000a000a0000000000000000,
            0x00000000000100010000000000000000,
            0x00000000000000000000000000110000,
            0x00000000000000000000000000850000,
            0x00000000000000c00000000000000010,
            0x00000000000000400000000000000080,
            0x00000000000000080004000000000000,
            0x00000000000000030003000000000000,
            0x00000000000000000002200000011000
        ],
        # Solution: [#6 of weight 135 and active sbox number 4]
        [
            0x00000000000000000000000000000a00,
            0x00000000000000000000000000000100,
            0x00000000000000010000000000000000,
            0x00000000000000080000000000000000,
            0x00080000000000000000000000000000,
            0x00030000000000000000000000000000,
            0x00000000000000002000000010000000,
            0x00000000000000006000000060000000,
            0x00004040000020200000000000000000,
            0x00005050000050500000000000000000,
            0x05050000000000000505000000000000,
            0x08020000000000000802000000000000,
            0x0000000000000000a000a00000000000,
            0x00000000000000001000100000000000,
            0x00000000000000000000110000000000,
            0x00000000000000000000880000000000,
            0x00000000000000000000080000000800,
            0x00000000000000000000030000000300,
            0x00000202000001010000000000000000,
            0x00000505000005050000000000000000,
            0x00000000050500000000000005050000,
            0x00000000080200000000000008020000,
            0x000000000000000000a000a000000000,
            0x00000000000000000010001000000000,
            0x00001100000000000000000000000000,
            0x00009800000000000000000000000000,
            0x00000000000000000900000008000000,
            0x00000000000000000300000003000000,
            0x00002020000010100000000000000000,
            0x00009060000080600000000000000000,
            0x0000000000000000050400000a0a0000,
            0x00000000000000000505000001010000,
            0x00000000000050100000000000005010,
            0x00000000000020800000000000002080,
            0x00000000000a000a0000000000000000,
            0x00000000000100010000000000000000,
            0x00000000000000000000000000110000,
            0x00000000000000000000000000850000,
            0x00000000000000c00000000000000010,
            0x00000000000000400000000000000080,
            0x00000000000000080004000000000000,
            0x00000000000000030003000000000000,
            0x00000000000000000002200000011000
        ],
        # No route with weight equal to 135.
    ],

],
