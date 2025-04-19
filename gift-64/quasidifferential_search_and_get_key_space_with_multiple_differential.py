import numpy, pyboolector

GIFT_SBOX = [0x1, 0xa, 0x4, 0xc, 0x6, 0xf, 0x3, 0x9, 0x2, 0xd, 0xb, 0x7, 0x5, 0x0, 0x8, 0xe]
n = 4
m = 4
sbox = GIFT_SBOX
state_bits = 64
state_words = 16
sbox_bits = 4

permutation_bits_table_64 = [
    0, 5, 10, 15, 16, 21, 26, 31, 32, 37, 42, 47, 48, 53, 58, 63,
    12, 1, 6, 11, 28, 17, 22, 27, 44, 33, 38, 43, 60, 49, 54, 59,
    8, 13, 2, 7, 24, 29, 18, 23, 40, 45, 34, 39, 56, 61, 50, 55,
    4, 9, 14, 3, 20, 25, 30, 19, 36, 41, 46, 35, 52, 57, 62, 51
]

round_constants = [0x01, 0x03, 0x07, 0x0F, 0x1F, 0x3E, 0x3D, 0x3B, 0x37, 0x2F, 0x1E, 0x3C, 0x39, 0x33, 0x27, 0x0E, 0x1D, 0x3A, 0x35, 0x2B, 0x16, 0x2C, 0x18, 0x30, 0x21, 0x02, 0x05, 0x0B, 0x17, 0x2E, 0x1C, 0x38]
round_constant_positions = [3, 7, 11, 15, 19, 23]

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

    for u in range(2 ** n):
        for v in range(2 ** m):
            count_x = 0
            for x in IN_a_to_b:
                c = vector_inner_product(u, v, x, sbox[x], n, m)
                count_x += (-1) ** c
            LAT[u, v] = count_x

    return LAT

def get_quasidifferentials_by_one_difference_route(difference_route, route_number, min_weight, max_weight):
    # min_weight = 0
    # max_weight = 1

    def search_quasidifferential(one_difference_characteristic, f):
        differential_characteristic_rounds = (len(one_difference_characteristic) - 1) // 2

        btor = pyboolector.Boolector()
        btor.Set_opt(pyboolector.BTOR_OPT_MODEL_GEN, 1)

        us = [btor.Var(btor.BitVecSort(state_bits), "u%d" % i) for i in range(differential_characteristic_rounds + 1)]
        vs = [btor.Var(btor.BitVecSort(state_bits), "v%d" % i) for i in range(differential_characteristic_rounds)]

        def permute_bits(x, y):
            for i in range(state_bits):
                btor.Assert(y[i] == x[permutation_bits_table_64[i]])

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
                        print("# error! a = {} b = {} has no routes".format(a, b), file=f)

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

                rcs = round_constants[round]
                mask_c = (1 << (state_bits - 1)) | (((rcs >> 5) & 0x1) << round_constant_positions[5]) | (
                            ((rcs >> 4) & 0x1) << round_constant_positions[4]) | (
                                     ((rcs >> 3) & 0x1) << round_constant_positions[3]) | (
                                     ((rcs >> 2) & 0x1) << round_constant_positions[2]) | (
                                     ((rcs >> 1) & 0x1) << round_constant_positions[1]) | (
                                     ((rcs >> 0) & 0x1) << round_constant_positions[0])
                if bin(mask_c & one_round_mask_s).count('1') % 2 == 1:
                    correlation_sign += 1

            return (-1) ** correlation_sign

        print("# differential : 0x{:016x} -> 0x{:016x}".format(one_difference_characteristic[0],
                                                               one_difference_characteristic[
                                                                   len(one_difference_characteristic) - 1]), file=f)
        quasidifferentials = []
        for target in range(min_weight, max_weight):
            mask_strings = []
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
                        print("     0x{:016x},".format(u), file=f)
                        print("     0x{:016x},".format(v), file=f)
                        trail.append(u)
                        trail.append(v)
                    s = compute_sign(one_difference_characteristic, trail)
                    if s < 0:
                        count_negative += 1
                    else:
                        count_positive += 1
                    trail.append(int(us[differential_characteristic_rounds].assignment, base=2))
                    mask_strings.append(trail)
                    print("     0x{:016x},".format(int(us[differential_characteristic_rounds].assignment, base=2)),
                          file=f)
                    print("     # Sign: {}".format(s), file=f)
                    print("    ],", file=f)
                    previous.append((s, [us[i].assignment for i in range(1, differential_characteristic_rounds)]))
                else:
                    print("     # No trails with weight equal to {}.".format(target), file=f)
                    print("     # sign = -1: {}".format(count_negative), file=f)
                    print("     # sign = +1: {}".format(count_positive), file=f)
                    break
            print("],", file=f)
            quasidifferentials.append(mask_strings)

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


import round_key_transition

def get_master_key_by_round_key(round, round_position, key_recovery_round_top):
    # key_recovery_round_top, key_recovery_round_bottom = get_key_recovery_extension_rounds()

    true_round = round + key_recovery_round_top

    if state_bits == 128:
        if round_position % 2 == 0:
            true_position = round_position // 2
            master_key = round_key_transition.key_recover_bit_from_top_to_bottom_128[true_round][true_position]
            return master_key
        elif round_position % 2 == 1:
            true_position = (round_position - 1) // 2
            master_key = round_key_transition.key_recover_bit_from_top_to_bottom_128[true_round][true_position]
            return master_key
    elif state_bits == 64:
        if round_position % 2 == 0:
            true_position = round_position // 2
            master_key = round_key_transition.key_recover_bit_from_top_to_bottom_64[true_round][true_position]
            return master_key
        elif round_position % 2 == 1:
            true_position = ((round_position - 1) // 2) + 1
            master_key = round_key_transition.key_recover_bit_from_top_to_bottom_64[true_round][true_position]
            return master_key


def get_key_space_by_one_difference_route_and_its_quasidifferential(difference_route, quasidifferentials, key_recovery_round_top):
    difference_characteristic_round = (len(difference_route) - 1) // 2

    def get_masterkey_matrix(difference_characteristic, linear_trail_all):

        one_trail_key_bit = 128
        masterkey_matrix1 = numpy.zeros((len(linear_trail_all), one_trail_key_bit + 1))
        masterkey_matrix = numpy.zeros((len(linear_trail_all), one_trail_key_bit))
        key_all_str = []
        count_positive_1 = 0
        count_negative_1 = 0
        for trail_number in range(len(linear_trail_all)):
            linear_trail = linear_trail_all[trail_number]

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
                                    if state_bits == 128:
                                        if k != 0 and k != 3:
                                            key_bit = (u >> k) & 0x1
                                            if key_bit == 1:
                                                round_key_full = round * state_bits + (sbox_bits * i + k)
                                                round_key_round = round_key_full // state_bits
                                                round_key_position = round_key_full % state_bits
                                                master_key = get_master_key_by_round_key(round_key_round,
                                                                                         round_key_position, key_recovery_round_top)
                                                master_key_word = master_key[0]
                                                master_key_position = master_key[1]
                                                master_key_full = master_key_word * 16 + master_key_position
                                                masterkey_matrix1[trail_number][master_key_full] = 1
                                                masterkey_matrix[trail_number][master_key_full] = 1
                                                correlation_key.append(
                                                    "k{}_{}".format(master_key_word, master_key_position))
                                    elif state_bits == 64:
                                        if k != 2 and k != 3:
                                            key_bit = (u >> k) & 0x1
                                            if key_bit == 1:
                                                round_key_full = round * state_bits + (sbox_bits * i + k)
                                                round_key_round = round_key_full // state_bits
                                                round_key_position = round_key_full % state_bits
                                                master_key = get_master_key_by_round_key(round_key_round,
                                                                                         round_key_position, key_recovery_round_top)
                                                master_key_word = master_key[0]
                                                master_key_position = master_key[1]
                                                master_key_full = master_key_word * 16 + master_key_position
                                                masterkey_matrix1[trail_number][master_key_full] = 1
                                                masterkey_matrix[trail_number][master_key_full] = 1
                                                correlation_key.append(
                                                    "k{}_{}".format(master_key_word, master_key_position))

                        c_weight = n - numpy.log2(abs(c))
                        if c < 0:
                            correlation_negative_1_number += 1

                        correlation *= c
                        correlation_weight += c_weight

                rcs = round_constants[round + key_recovery_round_top]
                mask_c = (1 << (state_bits - 1)) | (((rcs >> 5) & 0x1) << round_constant_positions[5]) | (
                            ((rcs >> 4) & 0x1) << round_constant_positions[4]) | (
                                     ((rcs >> 3) & 0x1) << round_constant_positions[3]) | (
                                     ((rcs >> 2) & 0x1) << round_constant_positions[2]) | (
                                     ((rcs >> 1) & 0x1) << round_constant_positions[1]) | (
                                     ((rcs >> 0) & 0x1) << round_constant_positions[0])
                correlation *= ((-1) ** bin(mask_c & one_round_mask_s).count('1'))
                if bin(mask_c & one_round_mask_s).count('1') % 2 == 1:
                    correlation_negative_1_number += 1

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
                one_trail_key_str = "linear equation: " + key_str
            else:
                one_trail_key_str = "linear equation: 1 = 1 (mask=0)"
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
        print("# ------------------------ expressions of key ------------------------")

        print('')
        print("# sign = -1 trail number = {}".format(count_negative_1))
        print("# sign = +1 trail number = {}".format(count_positive_1))

        return key_matrix, key_matrix1

    from sympy import Matrix
    def get_rank_and_base_master_key(key_matrix, key_matrix1):
        rank = numpy.linalg.matrix_rank(key_matrix)
        rank1 = numpy.linalg.matrix_rank(key_matrix1)
        print("# coefficient matrix's rank = {}, augmented matrix's rank = {}".format(rank, rank1))

        key_matrix_mat = Matrix(key_matrix1)
        print("# " + str(type(key_matrix_mat)))
        key_matrix_rref = numpy.array(key_matrix_mat.rref()[0].tolist()).astype(numpy.int32)
        key_rref_str = []
        key_space = []
        one_trail_key_bit = 128
        for r in range(rank1):
            one_trail_rref_str = ''
            one_key_linear_contraint = []
            for i in range(one_trail_key_bit):
                if key_matrix_rref[r][i] != 0:
                    word = i // 16
                    position = i % 16
                    one_trail_rref_str += "k" + str(word) + "[" + str(position) + "] + "
                    one_key_linear_contraint.append([word, position])
            one_trail_rref_str += " = " + str(key_matrix_rref[r][one_trail_key_bit])
            one_key_linear_contraint.append(key_matrix_rref[r][one_trail_key_bit])

            key_rref_str.append(one_trail_rref_str)
            key_space.append(one_key_linear_contraint)
        for k in range(len(key_rref_str)):
            print("# " + key_rref_str[k])

        return rank, rank1, key_rref_str, key_space

    key_matrix, key_matrix1 = print_all_trail_master_key()
    rank, rank1, master_key_space_str, master_key_space = get_rank_and_base_master_key(key_matrix, key_matrix1)

    return rank, rank1, master_key_space_str, master_key_space

from res import r_13_characteristics_p_64 as differentials
routes = differentials.routes[0]
print(len(routes))
min_weight = 0
max_weight = 1

all_valid_master_key_space = []
count_all_valid_master_key_space = []
count_all_valid_quasidifferential_number = []
count_all_valid_routes_index = []
count_all_valid_quasidifferential_number_with_routes_index = []
count_all_valid_differential_routes_index = []
count_all_valid_differential_quasidifferential_number = []
count_all_valid_differential_master_key_space = []

all_invalid_master_key_space = []
count_all_invalid_master_key_space = []
count_all_invalid_quasidifferential_number = []
count_all_invalid_routes_index = []

for rr in range(len(routes)):
    print("")
    difference_characteristic = routes[rr]
    quasidifferentials = get_quasidifferentials_by_one_difference_route(difference_characteristic, rr, min_weight, max_weight)
    print("route {} :".format(rr + 1))
    for w in range(len(quasidifferentials)):
        w_quasidifferentials = quasidifferentials[w]
        print("weight {} : {} quasidfferentials".format(w, len(w_quasidifferentials)))
        count_all_valid_quasidifferential_number_with_routes_index.append(len(w_quasidifferentials))
        rank, rank1, master_key_space_str, master_key_space = get_key_space_by_one_difference_route_and_its_quasidifferential(difference_characteristic, w_quasidifferentials, key_recovery_round_top=4)

        if rank == rank1:
            if master_key_space not in all_valid_master_key_space:
                all_valid_master_key_space.append(master_key_space)
                count_all_valid_master_key_space.append(1)
                count_all_valid_routes_index.append([rr + 1])
                count_all_valid_quasidifferential_number.append(len(w_quasidifferentials))

            elif master_key_space in all_valid_master_key_space:
                index_position = all_valid_master_key_space.index(master_key_space)
                count_all_valid_master_key_space[index_position] += 1
                routes_index_list = count_all_valid_routes_index[index_position]
                routes_index_list.append(rr + 1)
                count_all_valid_routes_index[index_position] = routes_index_list
                count_all_valid_quasidifferential_number[index_position] += len(w_quasidifferentials)
        elif rank != rank1:
            if master_key_space not in all_invalid_master_key_space:
                all_invalid_master_key_space.append(master_key_space)
                count_all_invalid_master_key_space.append(1)
                count_all_invalid_routes_index.append([rr + 1])
                count_all_invalid_quasidifferential_number.append(len(w_quasidifferentials))
            elif master_key_space in all_invalid_master_key_space:
                index_position = all_invalid_master_key_space.index(master_key_space)
                count_all_invalid_master_key_space[index_position] += 1
                routes_index_list = count_all_invalid_routes_index[index_position]
                routes_index_list.append(rr + 1)
                count_all_invalid_routes_index[index_position] = routes_index_list
                count_all_invalid_quasidifferential_number[index_position] += len(w_quasidifferentials)

print("")
print("{} valid possibilities of master key space".format(len(all_valid_master_key_space)))
for vp in range(len(all_valid_master_key_space)):
    print("")
    print("valid {} :".format(vp + 1))
    print("master key space : {}".format(all_valid_master_key_space[vp]))
    print("routes number in this space : {}".format(count_all_valid_master_key_space[vp]))
    print("routes index : {}".format(count_all_valid_routes_index[vp]))
    print("quasidifferential number in this space : {}".format(count_all_valid_quasidifferential_number[vp]))

    all_routes_index_in_this_master_key_space = count_all_valid_routes_index[vp]
    all_quasidifferential_number_with_routes_index_in_this_master_key_space = count_all_valid_quasidifferential_number_with_routes_index[vp]
    all_valid_output_difference_in_this_master_key_space = []
    all_valid_routes_index_with_output_difference_in_this_master_key_space = []
    all_valid_quasidifferential_number_with_output_difference_in_this_master_key_space = []
    for orr in range(len(all_routes_index_in_this_master_key_space)):
        one_route_index = all_routes_index_in_this_master_key_space[orr] - 1
        one_route = routes[one_route_index]
        quasidifferential_number_with_this_route = count_all_valid_quasidifferential_number_with_routes_index[one_route_index]

        one_differential = []
        one_differential.append(one_route[len(one_route) - 1])

        if one_differential not in all_valid_output_difference_in_this_master_key_space:
            all_valid_output_difference_in_this_master_key_space.append(one_differential)
            all_valid_routes_index_with_output_difference_in_this_master_key_space.append([one_route_index + 1])
            all_valid_quasidifferential_number_with_output_difference_in_this_master_key_space.append(quasidifferential_number_with_this_route)
        elif one_differential in all_valid_output_difference_in_this_master_key_space:
            output_difference_index_position = all_valid_output_difference_in_this_master_key_space.index(one_differential)
            routes_list = all_valid_routes_index_with_output_difference_in_this_master_key_space[output_difference_index_position]
            routes_list.append(one_route_index + 1)
            all_valid_routes_index_with_output_difference_in_this_master_key_space[output_difference_index_position] = routes_list
            all_valid_quasidifferential_number_with_output_difference_in_this_master_key_space[output_difference_index_position] += quasidifferential_number_with_this_route

    print("{} output differences in current space".format(len(all_valid_output_difference_in_this_master_key_space)))
    for dd in range(len((all_valid_output_difference_in_this_master_key_space))):
        print("      output difference {} : 0x{:016x}".format(dd + 1, all_valid_output_difference_in_this_master_key_space[dd][0]))
        print("      {}  routes index in this output difference : {}".format(len(all_valid_routes_index_with_output_difference_in_this_master_key_space[dd]), all_valid_routes_index_with_output_difference_in_this_master_key_space[dd]))
        print("      {} quasidifferential number in this output difference".format(all_valid_quasidifferential_number_with_output_difference_in_this_master_key_space[dd]))

print("")
print("{} invalid possibilities of master key space".format(len(all_invalid_master_key_space)))
for vp in range(len(all_invalid_master_key_space)):
    print("")
    print("invalid {} :".format(vp + 1))
    print("master key space : {}".format(all_invalid_master_key_space[vp]))
    print("routes number in this space : {}".format(count_all_invalid_master_key_space[vp]))
    print("routes index : {}".format(count_all_invalid_routes_index[vp]))
    print("quasidifferential number in this space : {}".format(count_all_invalid_quasidifferential_number[vp]))