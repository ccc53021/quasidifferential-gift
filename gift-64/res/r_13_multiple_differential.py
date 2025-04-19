# The 33 input differences of 33 characteristics, the same output difference is 0x0000000000001010,
# which are used to launched the 21-round multiple differential attack on GIFT-64.
# 16 with 4 quasidifferential trails satisfying |c| = p_{avg}
# 1 with 1 quasidifferential trails satisfying |c| = p_{avg}
# 16 with 16 quasidifferential trails satisfying |c| = p_{avg}
input_differences = [
    # 321
    # 3, 4
    0x000000c0000000c0,  # 46, 4
    0x000000c0000000d0,  # 53, 4
    0x000000c0000000f0,  # 57, 4
    0x000000e0000000c0,  # 62, 4
    0x000000d0000000c0,  # 75, 4
    0x000000c0000000e0,  # 85, 4
    0x000000f0000000c0,  # 86, 4
    0x000000d0000000e0,  # 116, 4
    0x000000e0000000d0,  # 131, 4
    0x000000e0000000f0,  # 136, 4
    0x000000f0000000d0,  # 144, 4
    0x000000d0000000d0,  # 146, 4
    0x000000e0000000e0,  # 150, 4
    0x000000f0000000e0,  # 155, 4
    0x000000d0000000f0,  # 180, 4
    0x000000f0000000f0,  # 186, 4
    # 4, 1
    0x000000f0000000d0,  # 4, 1
    # 7, 16
    0x000000f0000000f0,  # 41, 16
    0x000000c0000000c0,  # 50, 16
    0x000000c0000000f0,  # 51, 16
    0x000000c0000000d0,  # 54, 16
    0x000000d0000000f0,  # 58, 16
    0x000000f0000000e0,  # 59, 16
    0x000000c0000000e0,  # 60, 16
    0x000000e0000000c0,  # 102, 16
    0x000000e0000000f0,  # 110, 16
    0x000000e0000000d0,  # 111, 16
    0x000000d0000000e0,  # 112, 16
    0x000000d0000000c0,  # 162, 16
    0x000000f0000000d0,  # 163, 16
    0x000000d0000000d0,  # 181, 16
    0x000000f0000000c0,  # 191, 16
    0x000000e0000000e0,  # 192, 16
]

