import numpy, pyboolector

GIFT_SBOX = [0x1, 0xa, 0x4, 0xc, 0x6, 0xf, 0x3, 0x9, 0x2, 0xd, 0xb, 0x7, 0x5, 0x0, 0x8, 0xe]
n = 4
m = 4
sbox = GIFT_SBOX
state_bits = 128
state_words = 32
sbox_bits = 4

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

def get_DDT(sbox, n, m):
    DDT = numpy.zeros((2 ** n, 2 ** m))

    for input_d in range(2 ** n):
        for output_d in range(2 ** m):
            for x in range(2 ** n):
                y = x ^ input_d
                sx = sbox[x]
                sy = sbox[y]
                if (sx ^ sy) == output_d:
                    DDT[input_d, output_d] += 1

    return DDT


DDT = get_DDT(sbox, n, m)

def vector_inner_product(u, v, x, fx, n, m):
    left = 0
    for i in range(n):
        left += ((u >> i) & 0x1) * ((x >> i) & 0x1)
    left = left % 2
    right = 0
    for j in range(m):
        right += ((v >> j) & 0x1) * ((fx >> j) & 0x1)
    right = right % 2

    return left ^ right

def get_correlation_by_fixed_difference(a, b):
    IN_a_to_b = []
    for x in range(2 ** n):
        y = x ^ a
        sx = sbox[x]
        sy = sbox[y]
        if (sx ^ sy) == b:
            IN_a_to_b.append(x)

    LAT = numpy.zeros((2 ** n, 2 ** m))

    # (-1)^(<u,x>+<v,f(x)>)
    for u in range(2 ** n):
        for v in range(2 ** m):
            count_x = 0
            for x in IN_a_to_b:
                c = vector_inner_product(u, v, x, sbox[x], n, m)
                count_x += (-1) ** c
            LAT[u, v] = count_x

    return LAT


def get_key_recovery_extension_rounds(input_difference, output_difference):

    print(
        "# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
    top = distinguisher_extension.distinguisher_top_extension(input_difference)
    print("")
    print("# 0x{:032x} -> 0x{:032x}".format(input_difference, output_difference))
    print("")
    bottom = distinguisher_extension.distinguisher_bottom_extension(output_difference)

    return top, bottom


# quasidifferential one_solution_search
def get_quasidifferentials_by_one_difference_route(difference_route, route_number):
    min_weight = 0
    max_weight = 1

    def search_quasidifferential(one_difference_characteristic, f):
        differential_characteristic_rounds = (len(one_difference_characteristic) - 1) // 2

        btor = pyboolector.Boolector()
        btor.Set_opt(pyboolector.BTOR_OPT_MODEL_GEN, 1)

        # 轮掩码变量
        us = [btor.Var(btor.BitVecSort(state_bits), "u%d" % i) for i in range(differential_characteristic_rounds + 1)]
        vs = [btor.Var(btor.BitVecSort(state_bits), "v%d" % i) for i in range(differential_characteristic_rounds)]

        # 算法比特置换
        def permute_bits(x, y):
            for i in range(state_bits):
                btor.Assert(y[i] == x[permutation_bits_table_128[i]])

        def get_lat_weight(lat, flag):
            lat_weight = numpy.zeros((2 ** n, 2 ** m))
            if flag == 0:
                p = numpy.log2(abs(lat[0, 0]))
                for i in range(2 ** n):
                    for j in range(2 ** m):
                        c = lat[i, j]
                        if c != 0:
                            w = p - numpy.log2(abs(c))
                            lat_weight[i, j] = w
            else:
                if flag == 1:
                    if lat[0, 0] != 0:
                        p = numpy.log2(abs(lat[0, 0]))
                        for i in range(2 ** n):
                            for j in range(2 ** m):
                                c = lat[i, j]
                                if c != 0:
                                    w = p - numpy.log2(abs(c))
                                    lat_weight[i, j] = w
                                else:
                                    if c == 0:
                                        lat_weight[i, j] = 1
                    else:
                        print("# error! a = {} b = {} has no route".format(a, b), file=f)

            return lat_weight

        def get_one_words_weight(a, b, u, v):
            if a == b == 0:
                lat_weight = get_lat_weight(get_correlation_by_fixed_difference(0, 0), 0)
                weight0 = (u == 0) & (v == 0)

                weight1 = btor.Const(0)
                for x in range(2 ** n):
                    for y in range(2 ** m):
                        if lat_weight[x, y] == 1:
                            weight1 |= (u == x) & (v == y)

                weight2 = btor.Const(0)
                for x in range(2 ** n):
                    for y in range(2 ** m):
                        if lat_weight[x, y] == 2:
                            weight2 |= (u == x) & (v == y)

                btor.Assert(weight0 | weight1 | weight2)
                return btor.Cond(
                    weight1,
                    btor.Const(1, state_words),
                    btor.Cond(weight2, btor.Const(2, state_words), btor.Const(0, state_words))
                )
            else:
                lat_weight = get_lat_weight(get_correlation_by_fixed_difference(a, b), 1)
                allowed = btor.Const(0)
                for x in range(2 ** n):
                    for y in range(2 ** m):
                        if lat_weight[x, y] == 0:
                            allowed |= (u == x) & (v == y)
                btor.Assert(allowed)
                return btor.Const(0, state_words)

        cost = btor.Const(0, state_words)
        for i in range(differential_characteristic_rounds):
            permute_bits(vs[i], us[i + 1])
            for j in range(0, state_bits, 4):
                u = us[i][j + 3:j]
                v = vs[i][j + 3:j]
                a = (one_difference_characteristic[2 * i] >> j) & 0xf
                b = (one_difference_characteristic[2 * i + 1] >> j) & 0xf
                w = get_one_words_weight(a, b, u, v)
                cost += w

        btor.Assert(us[0] == 0)
        btor.Assert(us[differential_characteristic_rounds] == 0)
        btor.Set_opt(pyboolector.BTOR_OPT_INCREMENTAL, 1)

        def compute_sign(diffs, trail):
            correlation_sign = 0
            for round in range(differential_characteristic_rounds):
                one_round_difference_input = diffs[2 * round]
                one_round_difference_s = diffs[2 * round + 1]
                one_round_mask_input = trail[2 * round]
                one_round_mask_s = trail[2 * round + 1]

                for i in range(0, state_words):
                    a = (one_round_difference_input >> (4 * i)) & 0xf
                    b = (one_round_difference_s >> (4 * i)) & 0xf
                    u = (one_round_mask_input >> (4 * i)) & 0xf
                    v = (one_round_mask_s >> (4 * i)) & 0xf

                    if a != 0:
                        c = 1
                        if u == v == 0:
                            c = DDT[a, b]
                        else:
                            c = get_correlation_by_fixed_difference(a, b)[u, v]
                        c_weight = numpy.log2(abs(c)) - n
                        if c < 0:
                            correlation_sign += 1
            return (-1) ** correlation_sign

        print("# differential : 0x{:032x} -> 0x{:032x}".format(one_difference_characteristic[0],
                                                               one_difference_characteristic[
                                                                   len(one_difference_characteristic) - 1]), file=f)
        quasidifferentials = []
        for target in range(min_weight, max_weight):
            # mask_strings = []
            # Find all solutions
            previous = []
            print("# Solution: of weight {}".format(target), file=f)
            print("[", file=f)

            count_negative = 0
            count_positive = 0

            while True:
                btor.Assume(cost == target)
                distinct = btor.Const(1)
                for _, ws in previous:
                    temp = btor.Const(0)
                    for i in range(1, differential_characteristic_rounds):
                        temp |= (us[i] != btor.Const(ws[i - 1], state_bits))
                    distinct &= temp
                btor.Assume(distinct)

                r = btor.Sat()
                if r == btor.SAT:
                    print("    # Solution: [#{} of weight {}]".format(len(previous) + 1, target), file=f)
                    print("    [", file=f)
                    trail = []
                    for i in range(differential_characteristic_rounds):
                        u = int(us[i].assignment, base=2)
                        v = int(vs[i].assignment, base=2)
                        print("     0x{:032x},".format(u), file=f)
                        print("     0x{:032x},".format(v), file=f)
                        trail.append(u)
                        trail.append(v)
                    s = compute_sign(one_difference_characteristic, trail)
                    if s < 0:
                        count_negative += 1
                    else:
                        count_positive += 1
                    trail.append(int(us[differential_characteristic_rounds].assignment, base=2))
                    # mask_strings.append(trail)
                    quasidifferentials.append(trail)
                    print("     0x{:032x},".format(int(us[differential_characteristic_rounds].assignment, base=2)),
                          file=f)
                    print("     # Sign: {}".format(s), file=f)
                    print("    ],", file=f)
                    previous.append((s, [us[i].assignment for i in range(1, differential_characteristic_rounds)]))
                else:
                    print("     # No other trails with weight equal to {}.".format(target), file=f)
                    print("     # sign = -1: {}".format(count_negative), file=f)
                    print("     # sign = +1: {}".format(count_positive), file=f)
                    break
            print("],", file=f)
            # quasidifferentials.append(mask_strings)

        return quasidifferentials

    difference_characteristic_round = (len(difference_route) - 1) // 2
    filename = "result/quasidifferentials_r_{}_differential_route_{}.txt".format(difference_characteristic_round, route_number)
    f = open(filename, "a")
    print("# --------------------------------- route {} ---------------------------------------".format(route_number),
          file=f)
    print("routes = [".format(), file=f)
    quasidifferentials = search_quasidifferential(difference_route, f)
    print("],", file=f)
    print("", file=f)
    f.close()
    return quasidifferentials


import round_key_transition, distinguisher_extension

def get_master_key_by_round_key(round, position, key_recovery_round_top):

    # begin with round 1
    true_round = round + 1 + key_recovery_round_top
    true_position = position // 2
    u_or_v_position = true_position // 2
    full_position = u_or_v_position
    if true_position % 2 != 0:
        full_position += 32
    master_key = round_key_transition.master_key_in_round_from_top_to_bottom[true_round][full_position]

    return master_key


def get_key_space_by_one_difference_route(difference_route, quasidifferentials, key_recovery_round_top):
    difference_characteristic_round = (len(difference_route) - 1) // 2

    # get augmented matrix
    def get_masterkey_matrix(difference_characteristic, linear_trail_all):

        # one trail contains 128-bit master key
        one_trail_key_bit = 128
        # augmented matrix
        masterkey_matrix1 = numpy.zeros((len(linear_trail_all), one_trail_key_bit + 1))
        # coefficient
        masterkey_matrix = numpy.zeros((len(linear_trail_all), one_trail_key_bit))
        key_all_str = []
        count_positive_1 = 0
        count_negative_1 = 0
        for trail_number in range(len(linear_trail_all)):
            # one quasidifferential trail
            linear_trail = linear_trail_all[trail_number]

            # the correlation of the quasidifferential trail
            correlation = 1
            correlation_weight = 0
            correlation_negative_1_number = 0
            active_sbox = 0
            # linear expression of key
            correlation_key = []
            for round in range(0, difference_characteristic_round):

                one_round_difference_input = difference_characteristic[2 * round]
                one_round_difference_s = difference_characteristic[2 * round + 1]
                one_round_mask_input = linear_trail[2 * round]
                one_round_mask_s = linear_trail[2 * round + 1]

                for i in range(0, state_words):
                    a = (one_round_difference_input >> (sbox_bits * i)) & (2 ** sbox_bits - 1)
                    b = (one_round_difference_s >> (sbox_bits * i)) & (2 ** sbox_bits - 1)
                    u = (one_round_mask_input >> (sbox_bits * i)) & (2 ** sbox_bits - 1)
                    v = (one_round_mask_s >> (sbox_bits * i)) & (2 ** sbox_bits - 1)

                    if a != 0:
                        active_sbox += 1

                        if u == v == 0:
                            c = DDT[a, b]
                        else:
                            c = get_correlation_by_fixed_difference(a, b)[u, v]
                            if u != 0:
                                for k in range(sbox_bits):
                                    if k != 0 and k != 3:
                                        key_bit = (u >> k) & 0x1
                                        if key_bit == 1:
                                            # convert to master key
                                            round_key_full = round * state_bits + (sbox_bits * i + k)
                                            round_key_round = round_key_full // state_bits
                                            round_key_position = round_key_full % state_bits
                                            master_key = get_master_key_by_round_key(round_key_round,
                                                                                     round_key_position, key_recovery_round_top)
                                            master_key_round = master_key[0]
                                            master_key_position = master_key[1]
                                            master_key_full = master_key_round * 16 + master_key_position
                                            masterkey_matrix1[trail_number][master_key_full] = 1
                                            masterkey_matrix[trail_number][master_key_full] = 1
                                            correlation_key.append(
                                                "k{}_{}".format(master_key_round, master_key_position))
                        c_weight = n - numpy.log2(abs(c))
                        if c < 0:
                            correlation_negative_1_number += 1

                        correlation *= c
                        correlation_weight += c_weight

            if ((-1) ** correlation_negative_1_number) < 0:
                count_negative_1 += 1
            else:
                count_positive_1 += 1

            if len(correlation_key) > 0:
                key_str = ""
                for key in range(len(correlation_key) - 1):
                    key_str += correlation_key[key] + " + "
                key_str += correlation_key[len(correlation_key) - 1]
                if ((-1) ** correlation_negative_1_number) < 0:
                    key_str += " = 1"
                    masterkey_matrix1[trail_number][one_trail_key_bit] = 1
                else:
                    key_str += " = 0"
                one_trail_key_str = "linear equation：" + key_str
            else:
                one_trail_key_str = "linear equation：1 = 1 (mask=0)"
            one_line_str_len = 250
            for olsl in range(one_line_str_len - len(one_trail_key_str)):
                one_trail_key_str += " "
            one_trail_key_str += "(" + str(trail_number + 1) + ")"
            key_all_str.append(one_trail_key_str)

        return key_all_str, masterkey_matrix, masterkey_matrix1, count_positive_1, count_negative_1

    def print_all_trail_master_key():
        key_all_str, key_matrix, key_matrix1, count_positive_1, count_negative_1 = get_masterkey_matrix(
            difference_route, quasidifferentials)
        print('')
        print("# ------------------------ linear equations of key ------------------------")

        print('')
        print("# sign = -1 trail's number = {}".format(count_negative_1))
        print("# sign = +1 trail's number = {}".format(count_positive_1))

        return key_matrix, key_matrix1

    from sympy import Matrix

    def get_rank_and_base_master_key(key_matrix, key_matrix1):
        rank = numpy.linalg.matrix_rank(key_matrix)
        rank1 = numpy.linalg.matrix_rank(key_matrix1)
        print("# coefficient matrix's rank = {}, augmented matrix's rank = {}".format(rank, rank1))

        key_matrix_mat = Matrix(key_matrix1)
        key_matrix_rref = numpy.array(key_matrix_mat.rref()[0].tolist()).astype(numpy.int32)
        key_rref_str = []
        basis = []
        one_trail_key_bit = 128
        for r in range(rank1):
            one_b = []
            one_trail_rref_str = ''
            for i in range(one_trail_key_bit):
                if key_matrix_rref[r][i] != 0:
                    round = i // 16
                    position = i % 16
                    one_trail_rref_str += "k" + str(round) + "[" + str(position) + "] + "
                    one_b.append([round, position])
            one_trail_rref_str += " = " + str(key_matrix_rref[r][one_trail_key_bit])
            one_b.append([key_matrix_rref[r][one_trail_key_bit]])
            key_rref_str.append(one_trail_rref_str)
            basis.append(one_b)
        for k in range(len(key_rref_str)):
            print("# " + key_rref_str[k])

        return rank, rank1, key_rref_str, basis

    key_matrix, key_matrix1 = print_all_trail_master_key()
    rank, rank1, master_key_space, basis = get_rank_and_base_master_key(key_matrix, key_matrix1)

    return rank, rank1, master_key_space, basis

def get_key_recovery_extension_and_master_base_vector(route, qts):
    print("")
    print("")

    input_difference = route[0]
    output_difference = route[len(route) - 1]

    key_recovery_round_top, key_recovery_round_bottom = get_key_recovery_extension_rounds(input_difference,
                                                                                          output_difference)
    one_route_rank, one_route_rank1, one_route_master_key_space, basis = get_key_space_by_one_difference_route(route, qts, key_recovery_round_top)
    pro = numpy.log2(len(qts))
    print(one_route_master_key_space)
    print(basis)

    return key_recovery_round_top, one_route_rank, one_route_rank1, basis


def get_one_route_probability_by_base_vector(difference_route, quasidifferentials, key_recovery_round_top, base_vector):
    difference_characteristic_round = (len(difference_route) - 1) // 2

    def get_base_matrix_master_key(base_vector_num):
        matrix = numpy.zeros((len(base_vector_num), 129))
        for i in range(len(base_vector_num)):
            equation = base_vector_num[i]
            for e in range(len(equation) - 1):
                round = equation[e][0]
                position = equation[e][1]
                matrix[i][round * 16 + position] = 1
            matrix[i][128] = equation[len(equation) - 1][0]

        return matrix

    def get_one_trail_masterkey_vector(difference_characteristic, linear_trail):

        one_trail_key_bit = 128
        key_vector = [0] * (one_trail_key_bit + 1)

        correlation = 1
        correlation_weight = 0
        correlation_negative_1_number = 0
        active_sbox = 0
        correlation_key = []
        for round in range(0, difference_characteristic_round):

            one_round_difference_input = difference_characteristic[2 * round]
            one_round_difference_s = difference_characteristic[2 * round + 1]
            one_round_mask_input = linear_trail[2 * round]
            one_round_mask_s = linear_trail[2 * round + 1]

            for i in range(0, state_words):
                a = (one_round_difference_input >> (sbox_bits * i)) & (2 ** sbox_bits - 1)
                b = (one_round_difference_s >> (sbox_bits * i)) & (2 ** sbox_bits - 1)
                u = (one_round_mask_input >> (sbox_bits * i)) & (2 ** sbox_bits - 1)
                v = (one_round_mask_s >> (sbox_bits * i)) & (2 ** sbox_bits - 1)

                if a != 0:
                    active_sbox += 1

                    if u == v == 0:
                        c = DDT[a, b]
                    else:
                        c = get_correlation_by_fixed_difference(a, b)[u, v]
                        if u != 0:
                            for k in range(sbox_bits):
                                if k != 0 and k != 3:
                                    key_bit = (u >> k) & 0x1
                                    if key_bit == 1:
                                        round_key_full = round * state_bits + (sbox_bits * i + k)
                                        round_key_round = round_key_full // state_bits
                                        round_key_position = round_key_full % state_bits
                                        master_key = get_master_key_by_round_key(
                                            round_key_round,
                                            round_key_position, key_recovery_round_top)
                                        master_key_round = master_key[0]
                                        master_key_position = master_key[1]
                                        master_key_full = master_key_round * 16 + master_key_position
                                        correlation_key.append("k{}_{}".format(master_key_round, master_key_position))
                                        key_vector[master_key_full] = 1
                    c_weight = n - numpy.log2(abs(c))
                    if c < 0:
                        correlation_negative_1_number += 1

                    correlation *= c
                    correlation_weight += c_weight

        if len(correlation_key) > 0:
            key_str = ""
            for key in range(len(correlation_key) - 1):
                key_str += correlation_key[key] + " + "
            key_str += correlation_key[len(correlation_key) - 1]
            if ((-1) ** correlation_negative_1_number) < 0:
                key_str += " = 1"
                key_vector[one_trail_key_bit] = 1
            else:
                key_str += " = 0"
            one_trail_key_str = "linear equation：" + key_str
        else:
            # print("K = null（mask=0）")
            one_trail_key_str = "linear equation：1 = 1 (mask=0)"

        return one_trail_key_str, key_vector, (-1) ** correlation_negative_1_number

    def juge_new_trail_master_key(matrix, i):
        one_trail = quasidifferentials[i]
        one_trail_key_str, new_trail_key, sign = get_one_trail_masterkey_vector(difference_route, one_trail)
        rank = numpy.linalg.matrix_rank(matrix)
        matrix = numpy.row_stack((matrix, new_trail_key))
        new_rank = numpy.linalg.matrix_rank(matrix)

        flag = 0
        if rank == new_rank:
            flag = 1
        else:
            flag = 0

        return flag, sign

    def juge_all_trail_master_key(base_vector):
        matrix = get_base_matrix_master_key(base_vector)
        count_trail_number = 0
        count_sign_positive = 0
        count_sign_negative = 0
        for i in range(len(quasidifferentials)):
            flag, sign = juge_new_trail_master_key(matrix, i)
            count_trail_number += flag
            if flag == 1:
                if sign < 0:
                    count_sign_negative += 1
                else:
                    count_sign_positive += 1
        print('')
        print("quasidifferential's number is ：{}".format(len(quasidifferentials)))
        print("valid trail's number is : {}".format(count_trail_number))
        print("sign = -1 : {}, sign = +1 : {}".format(count_sign_negative, count_sign_positive))
        return count_trail_number, count_sign_negative, count_sign_positive

    juge_all_trail_master_key(base_vector)

# Take the 20-round differential 2 as an example, the procedure of 21-round differential is similar.
from res import r_20_differential_2 as diffs

difference_characteristics = []
for i in range(len(diffs.differential_2_20r)):
    diffs_p = diffs.differential_2_20r[i]
    for rj in diffs_p:
        difference_characteristics.append(rj)
quasidifferential_trails = []
for i in range(len(difference_characteristics)):
    route_i = difference_characteristics[i]
    qts_i = get_quasidifferentials_by_one_difference_route(route_i, i)
    quasidifferential_trails.append(qts_i)

def get_one_differential_probability():
    flag_y = 0
    for i in range(2):
        e0 = difference_characteristics[i]
        qts0 = quasidifferential_trails[i]

        key_recovery_round_top, one_route_rank, one_route_rank1, basis = get_key_recovery_extension_and_master_base_vector(e0, qts0)
        if one_route_rank == one_route_rank1:
            flag_y = 1

            for j in range(len(difference_characteristics)):
                ej = difference_characteristics[j]
                qtsj = quasidifferential_trails[j]
                get_one_route_probability_by_base_vector(ej, qtsj, key_recovery_round_top, basis)
            break

    if flag_y == 0:
        print("no usable difference characteristic")

get_one_differential_probability()
