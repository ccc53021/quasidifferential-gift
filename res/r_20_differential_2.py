# differential 2 in Tosc 2022, Table 7
# input difference is 0x000000000000000000000000000000a0
# output difference is 0x00000000000000002000000210000001
# p = 2^{-121.83} : 2 with 2^{-124}, 4 with 2^{-125}, 2 with 2^{-126}
# "Towards Key-recovery-attack Friendly Distinguishers: Application to GIFT-128
#       authors: Rui Zong, Xiaoyang Dong, Huaifeng Chen, Yiyuan Luo, Si Wang, and Zheng Li"
differential_2_20r = [
    # Solution: of weight 124
    [
        # Solution: [#1 of weight 124]
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
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x00000000000000000085000000000000,
            0x000000000000c0000000000000001000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001,
        ],
        # Solution: [#2 of weight 124]
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
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x000000000000000000a5000000000000,
            0x000000000000c0000000000000003000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001,
        ],
        # No route with weight equal to 124.
    ],
    # Solution: of weight 125
    [
        # Solution: [#1 of weight 125]
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
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x00000000000000000085000000000000,
            0x000000000000c0000000000000001000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001,
        ],
        # Solution: [#2 of weight 125]
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
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x000000000000000000a5000000000000,
            0x000000000000c0000000000000003000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001,
        ],
        # S,olution: [#3 of weight 125]
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
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x000000000000000000a5000000000000,
            0x000000000000c0000000000000003000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001,
        ],
        # Solution: [#4 of weight 125]
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
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x00000000000000000085000000000000,
            0x000000000000c0000000000000001000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001,
        ],
        # No route with weight equal to 125.
    ],
    # Solution: of weight 126
    [
        # Solution: [#1 of weight 126]
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
            0x00000098000000000000000000000000,
            0x09000000080000000000000000000000,
            0x03000000030000000000000000000000,
            0x20200000101000000000000000000000,
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x00000000000000000085000000000000,
            0x000000000000c0000000000000001000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001,
        ],
        # Solution: [#2 of weight 126]
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
            0x00000098000000000000000000000000,
            0x09000000080000000000000000000000,
            0x03000000030000000000000000000000,
            0x20200000101000000000000000000000,
            0x50500000505000000000000000000000,
            0x50500000000000005050000000000000,
            0x20800000000000002080000000000000,
            0x00000000a000a0000000000000000000,
            0x00000000100010000000000000000000,
            0x00000000000000000011000000000000,
            0x000000000000000000a5000000000000,
            0x000000000000c0000000000000003000,
            0x00000000000040000000000000008000,
            0x00040000000000000000000000000008,
            0x00030000000000000000000000000003,
            0x00000000000000002000000210000001,
        ],
        # No route with weight equal to 126.
    ],
    # Solution: of weight 127
    [
        # No route with weight equal to 127.
    ]
]