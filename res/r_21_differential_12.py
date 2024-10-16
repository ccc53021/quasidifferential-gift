# !/user/bin/env python3
# -*- coding: utf-8 -*-

# differential 12 in this paper
# input difference is 0x00000000000000000000000000a00000
# output difference is 0x00000000000000000000000408000002
# p = 2^{-129.64} : 130 * 0 + 131 * 1 + 132 * 2 + 133 * 1 + 134 * 1 + 135 * 3 = 129.64
differential_12_21r = [
    # Solution: of weight 130
    [
        # No route with weight equal to 130.
    ],
    # Solution: of weight 131
    [
        # Solution: [#1 of weight 131 and active sbox number 3]
        [
            0x00000000000000000000000000a00000,
            0x00000000000000000000000000100000,
            0x00000010000000000000000000000000,
            0x00000080000000000000000000000000,
            0x00000000080000000000000000000000,
            0x00000000030000000000000000000000,
            0x00200000001000000000000000000000,
            0x00600000006000000000000000000000,
            0x00000000000000004040000020200000,
            0x00000000000000005050000050500000,
            0x00005050000000000000505000000000,
            0x00008020000000000000802000000000,
            0x0000000000000000000000000a000a00,
            0x00000000000000000000000001000100,
            0x00000000000000110000000000000000,
            0x00000000000000880000000000000000,
            0x00080000000800000000000000000000,
            0x00030000000300000000000000000000,
            0x00000000000000002020000010100000,
            0x00000000000000005050000050500000,
            0x00005050000000000000505000000000,
            0x00008020000000000000802000000000,
            0x0000000000000000000000000a000a00,
            0x00000000000000000000000001000100,
            0x00000000000000110000000000000000,
            0x00000000000000980000000000000000,
            0x00090000000800000000000000000000,
            0x00010000000300000000000000000000,
            0x00000000000000000020000010100000,
            0x00000000000000000060000080600000,
            0x000000000000000000004040000020a0,
            0x00000000000000000000505000005010,
            0x00000505000000000000050100000000,
            0x00000208000000000000020800000000,
            0x0a000a00000000000000000000000000,
            0x01000100000000000000000000000000,
            0x00000000110000000000000000000000,
            0x000000009c0000000000000000000000,
            0x00000000000000000090000000c00000,
            0x00000000000000000010000000200000,
            0x00001000000000000000000000000020,
            0x00008000000000000000000000000060,
            0x00000000000000000000000408000002
        ],
        # No route with weight equal to 131.
    ],
    # Solution: of weight 132
    [
        # Solution: [#1 of weight 132 and active sbox number 3]
        [
            0x00000000000000000000000000a00000,
            0x00000000000000000000000000100000,
            0x00000010000000000000000000000000,
            0x00000080000000000000000000000000,
            0x00000000080000000000000000000000,
            0x00000000030000000000000000000000,
            0x00200000001000000000000000000000,
            0x00600000006000000000000000000000,
            0x00000000000000004040000020200000,
            0x00000000000000005050000050500000,
            0x00005050000000000000505000000000,
            0x00008020000000000000802000000000,
            0x0000000000000000000000000a000a00,
            0x00000000000000000000000001000100,
            0x00000000000000110000000000000000,
            0x00000000000000880000000000000000,
            0x00080000000800000000000000000000,
            0x00030000000300000000000000000000,
            0x00000000000000002020000010100000,
            0x00000000000000005050000050500000,
            0x00005050000000000000505000000000,
            0x00008020000000000000802000000000,
            0x0000000000000000000000000a000a00,
            0x00000000000000000000000001000100,
            0x00000000000000110000000000000000,
            0x00000000000000980000000000000000,
            0x00090000000800000000000000000000,
            0x00010000000300000000000000000000,
            0x00000000000000000020000010100000,
            0x00000000000000000060000090600000,
            0x000000000000000000004050000020a0,
            0x00000000000000000000505000005010,
            0x00000505000000000000050100000000,
            0x00000208000000000000020800000000,
            0x0a000a00000000000000000000000000,
            0x01000100000000000000000000000000,
            0x00000000110000000000000000000000,
            0x000000009c0000000000000000000000,
            0x00000000000000000090000000c00000,
            0x00000000000000000010000000200000,
            0x00001000000000000000000000000020,
            0x00008000000000000000000000000060,
            0x00000000000000000000000408000002
        ],
        # Solution: [#2 of weight 132 and active sbox number 3]
        [
            0x00000000000000000000000000a00000,
            0x00000000000000000000000000100000,
            0x00000010000000000000000000000000,
            0x00000080000000000000000000000000,
            0x00000000080000000000000000000000,
            0x00000000030000000000000000000000,
            0x00200000001000000000000000000000,
            0x00600000006000000000000000000000,
            0x00000000000000004040000020200000,
            0x00000000000000005050000050500000,
            0x00005050000000000000505000000000,
            0x00008020000000000000802000000000,
            0x0000000000000000000000000a000a00,
            0x00000000000000000000000001000100,
            0x00000000000000110000000000000000,
            0x00000000000000980000000000000000,
            0x00090000000800000000000000000000,
            0x00030000000300000000000000000000,
            0x00000000000000002020000010100000,
            0x00000000000000005050000050500000,
            0x00005050000000000000505000000000,
            0x00008020000000000000802000000000,
            0x0000000000000000000000000a000a00,
            0x00000000000000000000000001000100,
            0x00000000000000110000000000000000,
            0x00000000000000980000000000000000,
            0x00090000000800000000000000000000,
            0x00010000000300000000000000000000,
            0x00000000000000000020000010100000,
            0x00000000000000000060000080600000,
            0x000000000000000000004040000020a0,
            0x00000000000000000000505000005010,
            0x00000505000000000000050100000000,
            0x00000208000000000000020800000000,
            0x0a000a00000000000000000000000000,
            0x01000100000000000000000000000000,
            0x00000000110000000000000000000000,
            0x000000009c0000000000000000000000,
            0x00000000000000000090000000c00000,
            0x00000000000000000010000000200000,
            0x00001000000000000000000000000020,
            0x00008000000000000000000000000060,
            0x00000000000000000000000408000002
        ],
        # No route with weight equal to 132.
    ],
    # Solution: of weight 133
    [
        # Solution: [#1 of weight 133 and active sbox number 3]
        [
            0x00000000000000000000000000a00000,
            0x00000000000000000000000000100000,
            0x00000010000000000000000000000000,
            0x00000080000000000000000000000000,
            0x00000000080000000000000000000000,
            0x00000000030000000000000000000000,
            0x00200000001000000000000000000000,
            0x00600000006000000000000000000000,
            0x00000000000000004040000020200000,
            0x00000000000000005050000050500000,
            0x00005050000000000000505000000000,
            0x00008020000000000000802000000000,
            0x0000000000000000000000000a000a00,
            0x00000000000000000000000001000100,
            0x00000000000000110000000000000000,
            0x00000000000000980000000000000000,
            0x00090000000800000000000000000000,
            0x00030000000300000000000000000000,
            0x00000000000000002020000010100000,
            0x00000000000000005050000050500000,
            0x00005050000000000000505000000000,
            0x00008020000000000000802000000000,
            0x0000000000000000000000000a000a00,
            0x00000000000000000000000001000100,
            0x00000000000000110000000000000000,
            0x00000000000000980000000000000000,
            0x00090000000800000000000000000000,
            0x00010000000300000000000000000000,
            0x00000000000000000020000010100000,
            0x00000000000000000060000090600000,
            0x000000000000000000004050000020a0,
            0x00000000000000000000505000005010,
            0x00000505000000000000050100000000,
            0x00000208000000000000020800000000,
            0x0a000a00000000000000000000000000,
            0x01000100000000000000000000000000,
            0x00000000110000000000000000000000,
            0x000000009c0000000000000000000000,
            0x00000000000000000090000000c00000,
            0x00000000000000000010000000200000,
            0x00001000000000000000000000000020,
            0x00008000000000000000000000000060,
            0x00000000000000000000000408000002
        ],
        # No route with weight equal to 133.
    ],
    # Solution: of weight 134
    [
        # Solution: [#1 of weight 134 and active sbox number 3]
        [
            0x00000000000000000000000000a00000,
            0x00000000000000000000000000100000,
            0x00000010000000000000000000000000,
            0x00000080000000000000000000000000,
            0x00000000080000000000000000000000,
            0x00000000030000000000000000000000,
            0x00200000001000000000000000000000,
            0x00600000006000000000000000000000,
            0x00000000000000004040000020200000,
            0x00000000000000005050000050500000,
            0x00005050000000000000505000000000,
            0x00008020000000000000802000000000,
            0x0000000000000000000000000a000a00,
            0x00000000000000000000000001000100,
            0x00000000000000110000000000000000,
            0x00000000000000880000000000000000,
            0x00080000000800000000000000000000,
            0x00030000000300000000000000000000,
            0x00000000000000002020000010100000,
            0x00000000000000005050000050500000,
            0x00005050000000000000505000000000,
            0x00008020000000000000802000000000,
            0x0000000000000000000000000a000a00,
            0x00000000000000000000000001000100,
            0x00000000000000110000000000000000,
            0x00000000000000880000000000000000,
            0x00080000000800000000000000000000,
            0x00030000000300000000000000000000,
            0x00000000000000002020000010100000,
            0x00000000000000009060000080600000,
            0x0000000000000000000050400000a0a0,
            0x00000000000000000000505000001010,
            0x00000501000000000000050100000000,
            0x00000208000000000000020800000000,
            0x0a000a00000000000000000000000000,
            0x01000100000000000000000000000000,
            0x00000000110000000000000000000000,
            0x000000009c0000000000000000000000,
            0x00000000000000000090000000c00000,
            0x00000000000000000010000000200000,
            0x00001000000000000000000000000020,
            0x00008000000000000000000000000060,
            0x00000000000000000000000408000002
        ],
        # No route with weight equal to 134.
    ],
    # Solution: of weight 135
    [
        # Solution: [#1 of weight 135 and active sbox number 3]
        [
            0x00000000000000000000000000a00000,
            0x00000000000000000000000000100000,
            0x00000010000000000000000000000000,
            0x00000080000000000000000000000000,
            0x00000000080000000000000000000000,
            0x00000000030000000000000000000000,
            0x00200000001000000000000000000000,
            0x00600000006000000000000000000000,
            0x00000000000000004040000020200000,
            0x00000000000000005050000050500000,
            0x00005050000000000000505000000000,
            0x00008020000000000000802000000000,
            0x0000000000000000000000000a000a00,
            0x00000000000000000000000001000100,
            0x00000000000000110000000000000000,
            0x00000000000000880000000000000000,
            0x00080000000800000000000000000000,
            0x00030000000300000000000000000000,
            0x00000000000000002020000010100000,
            0x00000000000000005050000050500000,
            0x00005050000000000000505000000000,
            0x00008020000000000000802000000000,
            0x0000000000000000000000000a000a00,
            0x00000000000000000000000001000100,
            0x00000000000000110000000000000000,
            0x00000000000000880000000000000000,
            0x00080000000800000000000000000000,
            0x00030000000300000000000000000000,
            0x00000000000000002020000010100000,
            0x00000000000000009060000090600000,
            0x0000000000000000000050500000a0a0,
            0x00000000000000000000505000001010,
            0x00000501000000000000050100000000,
            0x00000208000000000000020800000000,
            0x0a000a00000000000000000000000000,
            0x01000100000000000000000000000000,
            0x00000000110000000000000000000000,
            0x000000009c0000000000000000000000,
            0x00000000000000000090000000c00000,
            0x00000000000000000010000000200000,
            0x00001000000000000000000000000020,
            0x00008000000000000000000000000060,
            0x00000000000000000000000408000002
        ],
        # Solution: [#2 of weight 135 and active sbox number 3]
        [
            0x00000000000000000000000000a00000,
            0x00000000000000000000000000100000,
            0x00000010000000000000000000000000,
            0x00000080000000000000000000000000,
            0x00000000080000000000000000000000,
            0x00000000030000000000000000000000,
            0x00200000001000000000000000000000,
            0x00600000006000000000000000000000,
            0x00000000000000004040000020200000,
            0x00000000000000005050000050500000,
            0x00005050000000000000505000000000,
            0x00008020000000000000802000000000,
            0x0000000000000000000000000a000a00,
            0x00000000000000000000000001000100,
            0x00000000000000110000000000000000,
            0x00000000000000880000000000000000,
            0x00080000000800000000000000000000,
            0x00030000000300000000000000000000,
            0x00000000000000002020000010100000,
            0x00000000000000005050000050500000,
            0x00005050000000000000505000000000,
            0x00008020000000000000802000000000,
            0x0000000000000000000000000a000a00,
            0x00000000000000000000000001000100,
            0x00000000000000110000000000000000,
            0x00000000000000980000000000000000,
            0x00090000000800000000000000000000,
            0x00030000000300000000000000000000,
            0x00000000000000002020000010100000,
            0x00000000000000009060000080600000,
            0x0000000000000000000050400000a0a0,
            0x00000000000000000000505000001010,
            0x00000501000000000000050100000000,
            0x00000208000000000000020800000000,
            0x0a000a00000000000000000000000000,
            0x01000100000000000000000000000000,
            0x00000000110000000000000000000000,
            0x000000009c0000000000000000000000,
            0x00000000000000000090000000c00000,
            0x00000000000000000010000000200000,
            0x00001000000000000000000000000020,
            0x00008000000000000000000000000060,
            0x00000000000000000000000408000002
        ],
        # Solution: [#3 of weight 135 and active sbox number 3]
        [
            0x00000000000000000000000000a00000,
            0x00000000000000000000000000100000,
            0x00000010000000000000000000000000,
            0x00000080000000000000000000000000,
            0x00000000080000000000000000000000,
            0x00000000030000000000000000000000,
            0x00200000001000000000000000000000,
            0x00600000006000000000000000000000,
            0x00000000000000004040000020200000,
            0x00000000000000005050000050500000,
            0x00005050000000000000505000000000,
            0x00008020000000000000802000000000,
            0x0000000000000000000000000a000a00,
            0x00000000000000000000000001000100,
            0x00000000000000110000000000000000,
            0x00000000000000980000000000000000,
            0x00090000000800000000000000000000,
            0x00030000000300000000000000000000,
            0x00000000000000002020000010100000,
            0x00000000000000005050000050500000,
            0x00005050000000000000505000000000,
            0x00008020000000000000802000000000,
            0x0000000000000000000000000a000a00,
            0x00000000000000000000000001000100,
            0x00000000000000110000000000000000,
            0x00000000000000880000000000000000,
            0x00080000000800000000000000000000,
            0x00030000000300000000000000000000,
            0x00000000000000002020000010100000,
            0x00000000000000009060000080600000,
            0x0000000000000000000050400000a0a0,
            0x00000000000000000000505000001010,
            0x00000501000000000000050100000000,
            0x00000208000000000000020800000000,
            0x0a000a00000000000000000000000000,
            0x01000100000000000000000000000000,
            0x00000000110000000000000000000000,
            0x000000009c0000000000000000000000,
            0x00000000000000000090000000c00000,
            0x00000000000000000010000000200000,
            0x00001000000000000000000000000020,
            0x00008000000000000000000000000060,
            0x00000000000000000000000408000002
        ],
        # No route with weight equal to 135.
    ],

],
