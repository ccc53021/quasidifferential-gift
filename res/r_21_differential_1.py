# !/user/bin/env python3
# -*- coding: utf-8 -*-

# differential 1 in this paper
# input difference is 0x000000000000000000000000000000a0
# output difference is 0x00000000000000002000000210000001
# p = 2^{-128.64} : 131 * 2 + 132 * 4 + 133 * 2 + 134 * 2 + 135 * 6 = 128.64
differential_1_21r = [
    # Solution: of weight 130
    [
        # No route with weight equal to 130.
    ],
    # Solution: of weight 131
    [
        # Solution: [#1 of weight 131 and active sbox number 4]
        [
            0x000000000000000000000000000000a0,
            0x00000000000000000000000000000010,
            0x00000001000000000000000000000000,
            0x00000008000000000000000000000000,
            0x08000000000000000000000000000000,
            0x03000000000000000000000000000000,
            0x20000000100000000000000000000000,
            0x60000000600000000000000000000000,
            0x40400000202000000000000000000000,
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x80200000000000008020000000000000,
            0x000000000000000000000000a000a000,
            0x00000000000000000000000010001000,
            0x00000000000000000000001100000000,
            0x00000000000000000000008800000000,
            0x00000800000008000000000000000000,
            0x00000300000003000000000000000000,
            0x02020000010100000000000000000000,
            0x05050000050500000000000000000000,
            0x00000000505000000000000050500000,
            0x00000000802000000000000080200000,
            0x00000000000000000000000000a000a0,
            0x00000000000000000000000000100010,
            0x00000011000000000000000000000000,
            0x0000008c000000000000000000000000,
            0x080000000c0000000000000000000000,
            0x03000000020000000000000000000000,
            0x20200000100000000000000000000000,
            0x60900000600000000000000000000000,
            0x50400000a02000000000000000000000,
            0x50500000105000000000000000000000,
            0x50100000000000005050000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x000000000000000000a5000000000000,
            0x000000000000c0000000000000003000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001
        ],
        # Solution: [#2 of weight 131 and active sbox number 4]
        [
            0x000000000000000000000000000000a0,
            0x00000000000000000000000000000010,
            0x00000001000000000000000000000000,
            0x00000008000000000000000000000000,
            0x08000000000000000000000000000000,
            0x03000000000000000000000000000000,
            0x20000000100000000000000000000000,
            0x60000000600000000000000000000000,
            0x40400000202000000000000000000000,
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x80200000000000008020000000000000,
            0x000000000000000000000000a000a000,
            0x00000000000000000000000010001000,
            0x00000000000000000000001100000000,
            0x00000000000000000000008800000000,
            0x00000800000008000000000000000000,
            0x00000300000003000000000000000000,
            0x02020000010100000000000000000000,
            0x05050000050500000000000000000000,
            0x00000000505000000000000050500000,
            0x00000000802000000000000080200000,
            0x00000000000000000000000000a000a0,
            0x00000000000000000000000000100010,
            0x00000011000000000000000000000000,
            0x0000008c000000000000000000000000,
            0x080000000c0000000000000000000000,
            0x03000000020000000000000000000000,
            0x20200000100000000000000000000000,
            0x60900000600000000000000000000000,
            0x50400000a02000000000000000000000,
            0x50500000105000000000000000000000,
            0x50100000000000005050000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x00000000000000000085000000000000,
            0x000000000000c0000000000000001000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001
        ],
        # No route with weight equal to 131.
    ],
    # Solution: of weight 132
    [
        # Solution: [#1 of weight 132 and active sbox number 4]
        [
            0x000000000000000000000000000000a0,
            0x00000000000000000000000000000010,
            0x00000001000000000000000000000000,
            0x00000008000000000000000000000000,
            0x08000000000000000000000000000000,
            0x03000000000000000000000000000000,
            0x20000000100000000000000000000000,
            0x60000000600000000000000000000000,
            0x40400000202000000000000000000000,
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x80200000000000008020000000000000,
            0x000000000000000000000000a000a000,
            0x00000000000000000000000010001000,
            0x00000000000000000000001100000000,
            0x00000000000000000000008800000000,
            0x00000800000008000000000000000000,
            0x00000300000003000000000000000000,
            0x02020000010100000000000000000000,
            0x05050000050500000000000000000000,
            0x00000000505000000000000050500000,
            0x00000000802000000000000080200000,
            0x00000000000000000000000000a000a0,
            0x00000000000000000000000000100010,
            0x00000011000000000000000000000000,
            0x0000009c000000000000000000000000,
            0x090000000c0000000000000000000000,
            0x03000000020000000000000000000000,
            0x20200000100000000000000000000000,
            0x60900000600000000000000000000000,
            0x50400000a02000000000000000000000,
            0x50500000105000000000000000000000,
            0x50100000000000005050000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x000000000000000000a5000000000000,
            0x000000000000c0000000000000003000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001
        ],
        # Solution: [#2 of weight 132 and active sbox number 4]
        [
            0x000000000000000000000000000000a0,
            0x00000000000000000000000000000010,
            0x00000001000000000000000000000000,
            0x00000008000000000000000000000000,
            0x08000000000000000000000000000000,
            0x03000000000000000000000000000000,
            0x20000000100000000000000000000000,
            0x60000000600000000000000000000000,
            0x40400000202000000000000000000000,
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x80200000000000008020000000000000,
            0x000000000000000000000000a000a000,
            0x00000000000000000000000010001000,
            0x00000000000000000000001100000000,
            0x00000000000000000000008800000000,
            0x00000800000008000000000000000000,
            0x00000300000003000000000000000000,
            0x02020000010100000000000000000000,
            0x05050000050500000000000000000000,
            0x00000000505000000000000050500000,
            0x00000000802000000000000080200000,
            0x00000000000000000000000000a000a0,
            0x00000000000000000000000000100010,
            0x00000011000000000000000000000000,
            0x0000009c000000000000000000000000,
            0x090000000c0000000000000000000000,
            0x03000000020000000000000000000000,
            0x20200000100000000000000000000000,
            0x60900000600000000000000000000000,
            0x50400000a02000000000000000000000,
            0x50500000105000000000000000000000,
            0x50100000000000005050000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x00000000000000000085000000000000,
            0x000000000000c0000000000000001000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001
        ],
        # Solution: [#3 of weight 132 and active sbox number 4]
        [
            0x000000000000000000000000000000a0,
            0x00000000000000000000000000000010,
            0x00000001000000000000000000000000,
            0x00000008000000000000000000000000,
            0x08000000000000000000000000000000,
            0x03000000000000000000000000000000,
            0x20000000100000000000000000000000,
            0x60000000600000000000000000000000,
            0x40400000202000000000000000000000,
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x80200000000000008020000000000000,
            0x000000000000000000000000a000a000,
            0x00000000000000000000000010001000,
            0x00000000000000000000001100000000,
            0x00000000000000000000009800000000,
            0x00000900000008000000000000000000,
            0x00000300000003000000000000000000,
            0x02020000010100000000000000000000,
            0x05050000050500000000000000000000,
            0x00000000505000000000000050500000,
            0x00000000802000000000000080200000,
            0x00000000000000000000000000a000a0,
            0x00000000000000000000000000100010,
            0x00000011000000000000000000000000,
            0x0000008c000000000000000000000000,
            0x080000000c0000000000000000000000,
            0x03000000020000000000000000000000,
            0x20200000100000000000000000000000,
            0x60900000600000000000000000000000,
            0x50400000a02000000000000000000000,
            0x50500000105000000000000000000000,
            0x50100000000000005050000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x000000000000000000a5000000000000,
            0x000000000000c0000000000000003000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001
        ],
        # Solution: [#4 of weight 132 and active sbox number 4]
        [
            0x000000000000000000000000000000a0,
            0x00000000000000000000000000000010,
            0x00000001000000000000000000000000,
            0x00000008000000000000000000000000,
            0x08000000000000000000000000000000,
            0x03000000000000000000000000000000,
            0x20000000100000000000000000000000,
            0x60000000600000000000000000000000,
            0x40400000202000000000000000000000,
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x80200000000000008020000000000000,
            0x000000000000000000000000a000a000,
            0x00000000000000000000000010001000,
            0x00000000000000000000001100000000,
            0x00000000000000000000009800000000,
            0x00000900000008000000000000000000,
            0x00000300000003000000000000000000,
            0x02020000010100000000000000000000,
            0x05050000050500000000000000000000,
            0x00000000505000000000000050500000,
            0x00000000802000000000000080200000,
            0x00000000000000000000000000a000a0,
            0x00000000000000000000000000100010,
            0x00000011000000000000000000000000,
            0x0000008c000000000000000000000000,
            0x080000000c0000000000000000000000,
            0x03000000020000000000000000000000,
            0x20200000100000000000000000000000,
            0x60900000600000000000000000000000,
            0x50400000a02000000000000000000000,
            0x50500000105000000000000000000000,
            0x50100000000000005050000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x00000000000000000085000000000000,
            0x000000000000c0000000000000001000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001
        ],
        # No route with weight equal to 132.
    ],
    # Solution: of weight 133
    [
        # Solution: [#1 of weight 133 and active sbox number 4]
        [
            0x000000000000000000000000000000a0,
            0x00000000000000000000000000000010,
            0x00000001000000000000000000000000,
            0x00000008000000000000000000000000,
            0x08000000000000000000000000000000,
            0x03000000000000000000000000000000,
            0x20000000100000000000000000000000,
            0x60000000600000000000000000000000,
            0x40400000202000000000000000000000,
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x80200000000000008020000000000000,
            0x000000000000000000000000a000a000,
            0x00000000000000000000000010001000,
            0x00000000000000000000001100000000,
            0x00000000000000000000009800000000,
            0x00000900000008000000000000000000,
            0x00000300000003000000000000000000,
            0x02020000010100000000000000000000,
            0x05050000050500000000000000000000,
            0x00000000505000000000000050500000,
            0x00000000802000000000000080200000,
            0x00000000000000000000000000a000a0,
            0x00000000000000000000000000100010,
            0x00000011000000000000000000000000,
            0x0000009c000000000000000000000000,
            0x090000000c0000000000000000000000,
            0x03000000020000000000000000000000,
            0x20200000100000000000000000000000,
            0x60900000600000000000000000000000,
            0x50400000a02000000000000000000000,
            0x50500000105000000000000000000000,
            0x50100000000000005050000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x000000000000000000a5000000000000,
            0x000000000000c0000000000000003000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001
        ],
        # Solution: [#2 of weight 133 and active sbox number 4]
        [
            0x000000000000000000000000000000a0,
            0x00000000000000000000000000000010,
            0x00000001000000000000000000000000,
            0x00000008000000000000000000000000,
            0x08000000000000000000000000000000,
            0x03000000000000000000000000000000,
            0x20000000100000000000000000000000,
            0x60000000600000000000000000000000,
            0x40400000202000000000000000000000,
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x80200000000000008020000000000000,
            0x000000000000000000000000a000a000,
            0x00000000000000000000000010001000,
            0x00000000000000000000001100000000,
            0x00000000000000000000009800000000,
            0x00000900000008000000000000000000,
            0x00000300000003000000000000000000,
            0x02020000010100000000000000000000,
            0x05050000050500000000000000000000,
            0x00000000505000000000000050500000,
            0x00000000802000000000000080200000,
            0x00000000000000000000000000a000a0,
            0x00000000000000000000000000100010,
            0x00000011000000000000000000000000,
            0x0000009c000000000000000000000000,
            0x090000000c0000000000000000000000,
            0x03000000020000000000000000000000,
            0x20200000100000000000000000000000,
            0x60900000600000000000000000000000,
            0x50400000a02000000000000000000000,
            0x50500000105000000000000000000000,
            0x50100000000000005050000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x00000000000000000085000000000000,
            0x000000000000c0000000000000001000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001
        ],
        # No route with weight equal to 133.
    ],
    # Solution: of weight 134
    [
        # Solution: [#1 of weight 134 and active sbox number 4]
        [
            0x000000000000000000000000000000a0,
            0x00000000000000000000000000000010,
            0x00000001000000000000000000000000,
            0x00000008000000000000000000000000,
            0x08000000000000000000000000000000,
            0x03000000000000000000000000000000,
            0x20000000100000000000000000000000,
            0x60000000600000000000000000000000,
            0x40400000202000000000000000000000,
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x80200000000000008020000000000000,
            0x000000000000000000000000a000a000,
            0x00000000000000000000000010001000,
            0x00000000000000000000001100000000,
            0x00000000000000000000008800000000,
            0x00000800000008000000000000000000,
            0x00000300000003000000000000000000,
            0x02020000010100000000000000000000,
            0x05050000050500000000000000000000,
            0x00000000505000000000000050500000,
            0x00000000802000000000000080200000,
            0x00000000000000000000000000a000a0,
            0x00000000000000000000000000100010,
            0x00000011000000000000000000000000,
            0x00000088000000000000000000000000,
            0x08000000080000000000000000000000,
            0x03000000030000000000000000000000,
            0x20200000101000000000000000000000,
            0x60900000608000000000000000000000,
            0x50400000a0a000000000000000000000,
            0x50500000101000000000000000000000,
            0x50100000000000005010000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x000000000000000000a5000000000000,
            0x000000000000c0000000000000003000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001
        ],
        # Solution: [#2 of weight 134 and active sbox number 4]
        [
            0x000000000000000000000000000000a0,
            0x00000000000000000000000000000010,
            0x00000001000000000000000000000000,
            0x00000008000000000000000000000000,
            0x08000000000000000000000000000000,
            0x03000000000000000000000000000000,
            0x20000000100000000000000000000000,
            0x60000000600000000000000000000000,
            0x40400000202000000000000000000000,
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x80200000000000008020000000000000,
            0x000000000000000000000000a000a000,
            0x00000000000000000000000010001000,
            0x00000000000000000000001100000000,
            0x00000000000000000000008800000000,
            0x00000800000008000000000000000000,
            0x00000300000003000000000000000000,
            0x02020000010100000000000000000000,
            0x05050000050500000000000000000000,
            0x00000000505000000000000050500000,
            0x00000000802000000000000080200000,
            0x00000000000000000000000000a000a0,
            0x00000000000000000000000000100010,
            0x00000011000000000000000000000000,
            0x00000088000000000000000000000000,
            0x08000000080000000000000000000000,
            0x03000000030000000000000000000000,
            0x20200000101000000000000000000000,
            0x60900000608000000000000000000000,
            0x50400000a0a000000000000000000000,
            0x50500000101000000000000000000000,
            0x50100000000000005010000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x00000000000000000085000000000000,
            0x000000000000c0000000000000001000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001
        ],
        # No route with weight equal to 134.
    ],
    # Solution: of weight 135
    [
        # Solution: [#1 of weight 135 and active sbox number 4]
        [
            0x000000000000000000000000000000a0,
            0x00000000000000000000000000000010,
            0x00000001000000000000000000000000,
            0x00000008000000000000000000000000,
            0x08000000000000000000000000000000,
            0x03000000000000000000000000000000,
            0x20000000100000000000000000000000,
            0x60000000600000000000000000000000,
            0x40400000202000000000000000000000,
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x80200000000000008020000000000000,
            0x000000000000000000000000a000a000,
            0x00000000000000000000000010001000,
            0x00000000000000000000001100000000,
            0x00000000000000000000008800000000,
            0x00000800000008000000000000000000,
            0x00000300000003000000000000000000,
            0x02020000010100000000000000000000,
            0x05050000050500000000000000000000,
            0x00000000505000000000000050500000,
            0x00000000802000000000000080200000,
            0x00000000000000000000000000a000a0,
            0x00000000000000000000000000100010,
            0x00000011000000000000000000000000,
            0x00000088000000000000000000000000,
            0x08000000080000000000000000000000,
            0x03000000030000000000000000000000,
            0x20200000101000000000000000000000,
            0x60900000609000000000000000000000,
            0x50500000a0a000000000000000000000,
            0x50500000101000000000000000000000,
            0x50100000000000005010000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x00000000000000000085000000000000,
            0x000000000000c0000000000000001000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001
        ],
        # Solution: [#2 of weight 135 and active sbox number 4]
        [
            0x000000000000000000000000000000a0,
            0x00000000000000000000000000000010,
            0x00000001000000000000000000000000,
            0x00000008000000000000000000000000,
            0x08000000000000000000000000000000,
            0x03000000000000000000000000000000,
            0x20000000100000000000000000000000,
            0x60000000600000000000000000000000,
            0x40400000202000000000000000000000,
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x80200000000000008020000000000000,
            0x000000000000000000000000a000a000,
            0x00000000000000000000000010001000,
            0x00000000000000000000001100000000,
            0x00000000000000000000009800000000,
            0x00000900000008000000000000000000,
            0x00000300000003000000000000000000,
            0x02020000010100000000000000000000,
            0x05050000050500000000000000000000,
            0x00000000505000000000000050500000,
            0x00000000802000000000000080200000,
            0x00000000000000000000000000a000a0,
            0x00000000000000000000000000100010,
            0x00000011000000000000000000000000,
            0x00000088000000000000000000000000,
            0x08000000080000000000000000000000,
            0x03000000030000000000000000000000,
            0x20200000101000000000000000000000,
            0x60900000608000000000000000000000,
            0x50400000a0a000000000000000000000,
            0x50500000101000000000000000000000,
            0x50100000000000005010000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x00000000000000000085000000000000,
            0x000000000000c0000000000000001000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001
        ],
        # Solution: [#3 of weight 135 and active sbox number 4]
        [
            0x000000000000000000000000000000a0,
            0x00000000000000000000000000000010,
            0x00000001000000000000000000000000,
            0x00000008000000000000000000000000,
            0x08000000000000000000000000000000,
            0x03000000000000000000000000000000,
            0x20000000100000000000000000000000,
            0x60000000600000000000000000000000,
            0x40400000202000000000000000000000,
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x80200000000000008020000000000000,
            0x000000000000000000000000a000a000,
            0x00000000000000000000000010001000,
            0x00000000000000000000001100000000,
            0x00000000000000000000008800000000,
            0x00000800000008000000000000000000,
            0x00000300000003000000000000000000,
            0x02020000010100000000000000000000,
            0x05050000050500000000000000000000,
            0x00000000505000000000000050500000,
            0x00000000802000000000000080200000,
            0x00000000000000000000000000a000a0,
            0x00000000000000000000000000100010,
            0x00000011000000000000000000000000,
            0x00000088000000000000000000000000,
            0x08000000080000000000000000000000,
            0x03000000030000000000000000000000,
            0x20200000101000000000000000000000,
            0x60900000609000000000000000000000,
            0x50500000a0a000000000000000000000,
            0x50500000101000000000000000000000,
            0x50100000000000005010000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x000000000000000000a5000000000000,
            0x000000000000c0000000000000003000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001
        ],
        # Solution: [#4 of weight 135 and active sbox number 4]
        [
            0x000000000000000000000000000000a0,
            0x00000000000000000000000000000010,
            0x00000001000000000000000000000000,
            0x00000008000000000000000000000000,
            0x08000000000000000000000000000000,
            0x03000000000000000000000000000000,
            0x20000000100000000000000000000000,
            0x60000000600000000000000000000000,
            0x40400000202000000000000000000000,
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x80200000000000008020000000000000,
            0x000000000000000000000000a000a000,
            0x00000000000000000000000010001000,
            0x00000000000000000000001100000000,
            0x00000000000000000000009800000000,
            0x00000900000008000000000000000000,
            0x00000300000003000000000000000000,
            0x02020000010100000000000000000000,
            0x05050000050500000000000000000000,
            0x00000000505000000000000050500000,
            0x00000000802000000000000080200000,
            0x00000000000000000000000000a000a0,
            0x00000000000000000000000000100010,
            0x00000011000000000000000000000000,
            0x00000088000000000000000000000000,
            0x08000000080000000000000000000000,
            0x03000000030000000000000000000000,
            0x20200000101000000000000000000000,
            0x60900000608000000000000000000000,
            0x50400000a0a000000000000000000000,
            0x50500000101000000000000000000000,
            0x50100000000000005010000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x000000000000000000a5000000000000,
            0x000000000000c0000000000000003000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001
        ],
        # Solution: [#5 of weight 135 and active sbox number 4]
        [
            0x000000000000000000000000000000a0,
            0x00000000000000000000000000000010,
            0x00000001000000000000000000000000,
            0x00000008000000000000000000000000,
            0x08000000000000000000000000000000,
            0x03000000000000000000000000000000,
            0x20000000100000000000000000000000,
            0x60000000600000000000000000000000,
            0x40400000202000000000000000000000,
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x80200000000000008020000000000000,
            0x000000000000000000000000a000a000,
            0x00000000000000000000000010001000,
            0x00000000000000000000001100000000,
            0x00000000000000000000008800000000,
            0x00000800000008000000000000000000,
            0x00000300000003000000000000000000,
            0x02020000010100000000000000000000,
            0x05050000050500000000000000000000,
            0x00000000505000000000000050500000,
            0x00000000802000000000000080200000,
            0x00000000000000000000000000a000a0,
            0x00000000000000000000000000100010,
            0x00000011000000000000000000000000,
            0x00000098000000000000000000000000,
            0x09000000080000000000000000000000,
            0x03000000030000000000000000000000,
            0x20200000101000000000000000000000,
            0x60900000608000000000000000000000,
            0x50400000a0a000000000000000000000,
            0x50500000101000000000000000000000,
            0x50100000000000005010000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x00000000000000000085000000000000,
            0x000000000000c0000000000000001000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001
        ],
        # Solution: [#6 of weight 135 and active sbox number 4]
        [
            0x000000000000000000000000000000a0,
            0x00000000000000000000000000000010,
            0x00000001000000000000000000000000,
            0x00000008000000000000000000000000,
            0x08000000000000000000000000000000,
            0x03000000000000000000000000000000,
            0x20000000100000000000000000000000,
            0x60000000600000000000000000000000,
            0x40400000202000000000000000000000,
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x80200000000000008020000000000000,
            0x000000000000000000000000a000a000,
            0x00000000000000000000000010001000,
            0x00000000000000000000001100000000,
            0x00000000000000000000008800000000,
            0x00000800000008000000000000000000,
            0x00000300000003000000000000000000,
            0x02020000010100000000000000000000,
            0x05050000050500000000000000000000,
            0x00000000505000000000000050500000,
            0x00000000802000000000000080200000,
            0x00000000000000000000000000a000a0,
            0x00000000000000000000000000100010,
            0x00000011000000000000000000000000,
            0x00000098000000000000000000000000,
            0x09000000080000000000000000000000,
            0x03000000030000000000000000000000,
            0x20200000101000000000000000000000,
            0x60900000608000000000000000000000,
            0x50400000a0a000000000000000000000,
            0x50500000101000000000000000000000,
            0x50100000000000005010000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x000000000000000000a5000000000000,
            0x000000000000c0000000000000003000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001
        ],
        # No route with weight equal to 135.
    ],
]
