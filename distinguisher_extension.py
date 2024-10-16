# -*- coding: utf-8 -*-

state_bits = 128

permutation_bits_table_128 = [
    0, 5, 10, 15, 16, 21, 26, 31, 32, 37, 42, 47, 48, 53, 58, 63,
    64, 69, 74, 79, 80, 85, 90, 95, 96, 101, 106, 111, 112, 117, 122, 127,
    12, 1, 6, 11, 28, 17, 22, 27, 44, 33, 38, 43, 60, 49, 54, 59,
    76, 65, 70, 75, 92, 81, 86, 91, 108, 97, 102, 107, 124, 113, 118, 123,
    8, 13, 2, 7, 24, 29, 18, 23, 40, 45, 34, 39, 56, 61, 50, 55,
    72, 77, 66, 71, 88, 93, 82, 87, 104, 109, 98, 103, 120, 125, 114, 119,
    4, 9, 14, 3, 20, 25, 30, 19, 36, 41, 46, 35, 52, 57, 62, 51,
    68, 73, 78, 67, 84, 89, 94, 83, 100, 105, 110, 99, 116, 121, 126, 115
]
permutation_bits_table_128_inverse = [
    0, 33, 66, 99, 96, 1, 34, 67, 64, 97, 2, 35, 32, 65, 98, 3,
    4, 37, 70, 103, 100, 5, 38, 71, 68, 101, 6, 39, 36, 69, 102, 7,
    8, 41, 74, 107, 104, 9, 42, 75, 72, 105, 10, 43, 40, 73, 106, 11,
    12, 45, 78, 111, 108, 13, 46, 79, 76, 109, 14, 47, 44, 77, 110, 15,
    16, 49, 82, 115, 112, 17, 50, 83, 80, 113, 18, 51, 48, 81, 114, 19,
    20, 53, 86, 119, 116, 21, 54, 87, 84, 117, 22, 55, 52, 85, 118, 23,
    24, 57, 90, 123, 120, 25, 58, 91, 88, 121, 26, 59, 56, 89, 122, 27,
    28, 61, 94, 127, 124, 29, 62, 95, 92, 125, 30, 63, 60, 93, 126, 31
]

def get_difference_p_inverse(difference_1d):
    result = [0] * state_bits
    for i in range(state_bits):
        result[i] = difference_1d[permutation_bits_table_128_inverse[i]]

    return result

def get_difference_p(difference_1d):
    result = [0] * state_bits
    for i in range(state_bits):
        result[i] = difference_1d[permutation_bits_table_128[i]]

    return result

def get_difference_s_inverse(difference_1d):
    # 1*** -> 0001
    # *1** -> 0010
    # 11** -> 0100
    result = [0] * state_bits
    for i in range(0, state_bits, 4):
        if (difference_1d[i + 0] == 1) and (difference_1d[i + 1] == 0) and (difference_1d[i + 2] == 0) and (difference_1d[i + 3] == 0):
            result[i + 0] = 2
            result[i + 1] = 2
            result[i + 2] = 2
            result[i + 3] = 1
        elif (difference_1d[i + 0] == 0) and (difference_1d[i + 1] == 1) and (difference_1d[i + 2] == 0) and (difference_1d[i + 3] == 0):
            result[i + 0] = 2
            result[i + 1] = 2
            result[i + 2] = 1
            result[i + 3] = 2
        elif (difference_1d[i + 0] == 0) and (difference_1d[i + 1] == 0) and (difference_1d[i + 2] == 1) and (difference_1d[i + 3] == 0):
            result[i + 0] = 2
            result[i + 1] = 2
            result[i + 2] = 1
            result[i + 3] = 1
        elif (difference_1d[i + 0] == 0) and (difference_1d[i + 1] == 0) and (difference_1d[i + 2] == 0) and (difference_1d[i + 3] == 0):
            result[i + 0] = 0
            result[i + 1] = 0
            result[i + 2] = 0
            result[i + 3] = 0
        else:
            result[i + 0] = 2
            result[i + 1] = 2
            result[i + 2] = 2
            result[i + 3] = 2

    return result

def get_difference_s(difference_1d):
    # 0100 -> ***1
    # 1000 -> **11
    result = [0] * state_bits
    for i in range(0, state_bits, 4):
        if (difference_1d[i + 0] == 0) and (difference_1d[i + 1] == 0) and (difference_1d[i + 2] == 1) and (difference_1d[i + 3] == 0):
            result[i + 0] = 1
            result[i + 1] = 2
            result[i + 2] = 2
            result[i + 3] = 2
        elif (difference_1d[i + 0] == 0) and (difference_1d[i + 1] == 0) and (difference_1d[i + 2] == 0) and (difference_1d[i + 3] == 1):
            result[i + 0] = 1
            result[i + 1] = 1
            result[i + 2] = 2
            result[i + 3] = 2
        elif (difference_1d[i + 0] == 0) and (difference_1d[i + 1] == 0) and (difference_1d[i + 2] == 0) and (difference_1d[i + 3] == 0):
            result[i + 0] = 0
            result[i + 1] = 0
            result[i + 2] = 0
            result[i + 3] = 0
        else:
            result[i + 0] = 2
            result[i + 1] = 2
            result[i + 2] = 2
            result[i + 3] = 2

    return result

def judge_state(difference_states):
    count_2 = difference_states.count(2)
    # print(count_2)
    if count_2 == state_bits:
        return 1, count_2
    else:
        return 0, count_2

def get_states_by_difference(difference_1d):
    states = [0] * state_bits
    for i in range(state_bits):
        if difference_1d[i] == 0:
            states[i] = 0
        elif difference_1d[i] == 1:
            states[i] = 1
        else:
            states[i] = 2

    return states

def format_output(difference_1d):
    output_str = ""
    for i in range(0, state_bits, 4):
        for j in range(4):
            if difference_1d[i + j] == 0:
                output_str += "-"
            elif difference_1d[i + j] == 1:
                output_str += "1"
            elif difference_1d[i + j] == 2:
                output_str += "*"
        output_str += " "

    print("# " + output_str[::-1])

    return output_str

def distinguisher_top_extension(distinguisher_input_difference):
    # RK, P^-1, S^-1
    difference_initial = [0] * state_bits
    for i in range(state_bits):
        difference_initial[i] = (distinguisher_input_difference >> i) & 0x1
    difference_x_s = difference_initial

    flag_full = 0
    difference_list = []
    while (flag_full == 0):
        difference_x_p = get_difference_p_inverse(difference_x_s)
        # format_output(difference_x_p)
        states_x_p = get_states_by_difference(difference_x_p)
        p_list = [difference_x_p, states_x_p]

        difference_x_s = get_difference_s_inverse(difference_x_p)
        # format_output(difference_x_s)
        states_x_s = get_states_by_difference(difference_x_s)
        s_list = [difference_x_s, states_x_s]

        difference_list.append([p_list, s_list])

        flag_full, count_2 = judge_state(states_x_s)
        # print(flag_full)
    round = len(difference_list)
    print("# ------------------------ top: {} rounds extension -------------------------".format(round))
    for i in range(round):
        one_round = difference_list[round - i - 1]
        p_list = one_round[0]
        s_list = one_round[1]
        difference_x_p = p_list[0]
        difference_x_s = s_list[0]
        states_x_s = s_list[1]
        flag_full, count_2 = judge_state(states_x_s)
        print("# round {}: unknown {} bits".format(i + 1, count_2))
        format_output(difference_x_s)
        format_output(difference_x_p)

    return round

def distinguisher_bottom_extension(distinguisher_output_difference):
    # RK', S^-1, P^-1
    difference_initial = [0] * state_bits
    for i in range(state_bits):
        difference_initial[i] = (distinguisher_output_difference >> i) & 0x1
    difference_x_p = difference_initial

    flag_full = 0
    difference_list = []
    while (flag_full == 0):
        difference_x_s = get_difference_s(difference_x_p)
        # format_output(difference_x_s)
        states_x_s = get_states_by_difference(difference_x_s)
        s_list = [difference_x_s, states_x_s]

        difference_x_p = get_difference_p(difference_x_s)
        states_x_p = get_states_by_difference(difference_x_p)
        p_list = [difference_x_p, states_x_p]

        difference_list.append([s_list, p_list])

        flag_full, count_2 = judge_state(states_x_s)
        # print(flag_full)
    round = len(difference_list) - 1
    print("# ------------------------ bottom: {} rounds extension -------------------------".format(round))
    for i in range(round):
        one_round = difference_list[i]
        s_list = one_round[0]
        p_list = one_round[1]
        difference_x_s = s_list[0]
        difference_x_p = p_list[0]
        states_x_p = p_list[1]
        flag_full, count_2 = judge_state(states_x_p)
        print("# round {}: unknown {} bits".format(i + 1, count_2))
        format_output(difference_x_s)
        format_output(difference_x_p)

    return round


# distinguisher_top_extension(0x00000000000000000000000000000a00)
# distinguisher_bottom_extension(0x00140000008200000000000000000000)
