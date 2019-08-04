import random
import unittest

import yaml
from itertools import combinations

import linearcode
from linearcode import LinearCode, \
    add_matrix_in_front, \
    allocate_bits_in_column, \
    calculate_p_error, \
    compute_parity, \
    f_l_i_v_inner_loop, \
    fill_generator_matrix_from_a_matrix_transposed, \
    fill_i_matrix, \
    fill_parity_check_matrix_and_a_matrix_transposed, \
    find_linearly_independent_vector, \
    generate_syndrome_decoding_table, \
    get_binomial_coefficient, \
    get_hamming_weight, \
    is_gilbert_varshamov_bound, \
    is_power_of_two, \
    multiply_matrices, \
    partition_binary_numbers_by_weight, \
    partition_by_is_power_of_two, \
    read_file_to_dict, \
    sum_modulo_2, \
    transpose_matrix, \
    get_random_number_of_hamming_weight


# TODO test exceptions
# TODO test files
# TODO test ErrorCorrectingCode
class LinearCodeTest(unittest.TestCase):
    n = 21
    k = 6
    r = 15
    d = 7

    def test_init(self):
        LinearCode(n=self.n, k=self.k, d=self.d, channel_error_probability=0.01)
        pass

    def test_is_gilbert_varshamov_bound(self):
        #       /10 - 1\     /10 - 1\             5
        # 1  +  |      |  +  |      |  =  46  >  2  =  32  => false
        #       \   1  /     \   2  /
        assert not is_gilbert_varshamov_bound(10, 5, 4)
        #       /12 - 1\     /12 - 1\             8
        # 1  +  |      |  +  |      |  =  67  <  2  =  256  => true
        #       \   1  /     \   2  /
        assert is_gilbert_varshamov_bound(12, 4, 4)
        pass

    def test_find_linearly_independent_vector(self):
        assert is_gilbert_varshamov_bound(self.n, self.k, self.d)
        i_matrix = fill_i_matrix(self.r)
        for i in range(1000):
            result = find_linearly_independent_vector(i_matrix, self.r, self.d - 2)
            assert get_hamming_weight(result) >= self.d - 1
        pass

    def test_f_l_i_v_inner_loop(self):
        linear_combination_size = self.d - 2
        combinations_v_p = []
        vector_pool = fill_i_matrix(self.r)
        for i in range(linear_combination_size):
            combinations_v_p.append(list(combinations(vector_pool, i + 1)))
        test_value = 4378
        assert get_hamming_weight(test_value) < self.d - 1
        for i in range(1000):
            assert f_l_i_v_inner_loop(test_value, combinations_v_p) != test_value

    def test_sum_modulo_2(self):
        test_list = []
        n = 10
        print()
        for i in range(n):
            test_list.append(i + 1)
            print('{0:0>{width}b}'.format(i + 1, width=n))
        assert sum_modulo_2(iterable=test_list) == 11

    def test_allocate_bits_in_column(self):
        matrix = [0, 0, 0]
        n = 10
        position = 1
        allocate_bits_in_column(matrix=matrix, position=position, number=2, n=n)
        assert 2 ** (n - 1 - position) == matrix[1]
        pass

    def test_fill_parity_check_matrix_and_a_matrix_transposed(self):
        parity_check_matrix, a_matrix_transposed = fill_parity_check_matrix_and_a_matrix_transposed(
            self.n,
            self.r,
            self.d)
        print()
        for line in parity_check_matrix:
            print(format(line, '#0{0}b'.format(n + 2)))
        pass

    def test_is_power_of_two(self):
        num = 256
        assert is_power_of_two(num)
        num = 255
        assert not is_power_of_two(num)
        pass

    def test_fill_i_matrix(self):
        r = 10
        matrix = fill_i_matrix(size=r)
        i = 0
        for num in matrix:
            assert 2 ** (r - 1 - i) == num
            i += 1
        pass

    def test_append_matrix_to_front(self):
        n = 10
        k = 7
        r = n - k
        matrix1 = []
        print()
        for i in range(k):
            matrix1.append(i + 1)
        matrix2 = fill_i_matrix(k)
        i = 0
        for num in add_matrix_in_front(matrix1, matrix2, r):
            # each number consists of the power of two and an element of an I-matrix
            assert num == 2 ** (n - 1 - i) + (i + 1)
            i += 1
        pass

    def test_generator_matrix_from_a_matrix_transposed(self):
        n = 10
        k = 3
        r = n - k
        d = 1
        parity_check_matrix, a_matrix_transposed = fill_parity_check_matrix_and_a_matrix_transposed(n, r, d)
        print()
        for i in parity_check_matrix:
            print(format(i, '#0{0}b'.format(n + 2)))
        print()
        for i in fill_generator_matrix_from_a_matrix_transposed(a_matrix_transposed, k, r):
            print(format(i, '#0{0}b'.format(n + 2)))

    def test_transpose_matrix(self):
        n = 10
        r = 7
        d = 1
        matrix, a_matrix_transposed = fill_parity_check_matrix_and_a_matrix_transposed(n, r, d)
        print()
        for i in matrix:
            print(format(i, '#0{0}b'.format(n + 2)))
        print()
        transposed_matrix = transpose_matrix(matrix, n)
        for i in transposed_matrix:
            print(format(i, '#0{0}b'.format(r + 2)))
        print()
        doubly_transposed_matrix = transpose_matrix(transposed_matrix, r)

        for i in range(len(doubly_transposed_matrix)):
            assert doubly_transposed_matrix[i] == matrix[i]
            print(format(doubly_transposed_matrix[i], '#0{0}b'.format(n + 2)))

    def test_compute_parity(self):
        number_parity_0 = 0b101
        number_parity_exceeding_length_0 = 0b1101
        number_parity_1 = 0b11111111111
        number_parity_exceeding_length_1 = 0b111111111111
        assert compute_parity(number_parity_0, 3) == 0
        assert compute_parity(number_parity_exceeding_length_0, 3) == 0
        assert compute_parity(number_parity_1, 11) == 1
        assert compute_parity(number_parity_exceeding_length_1, 3) == 1
        pass

    def test_multiply_matrices(self):
        parity_check_matrix, a_matrix_transposed = fill_parity_check_matrix_and_a_matrix_transposed(
            self.n,
            self.r,
            self.d)
        generator_matrix = fill_generator_matrix_from_a_matrix_transposed(a_matrix_transposed, self.k, self.r)
        multiplied_matrix = multiply_matrices(
            parity_check_matrix,
            self.n,
            transpose_matrix(generator_matrix, self.n),
            self.k)
        print()
        for i in multiplied_matrix:
            assert i == 0
        pass

    def test_partition_by_is_power_of_two(self):
        n = 100
        powers_of_two, not_powers_of_two = partition_by_is_power_of_two(n)
        max_value = 0
        for i in range(n):
            max_value += i + 1
        i = 1
        counter = 0
        while True:
            if i <= n:
                counter += 1
                i *= 2
            else:
                break
        assert len(powers_of_two) == counter
        assert len(not_powers_of_two) == n - counter
        assert sum(powers_of_two) + sum(not_powers_of_two) == max_value
        pass

    def test_get_hamming_weight(self):
        assert get_hamming_weight(7) == 3
        assert get_hamming_weight(8) == 1
        assert get_hamming_weight(1) == 1
        pass

    def test_get_binary_numbers_partitioned_by_weight(self):
        n = 127
        length = len(bin(n)) - 2
        binary_numbers_partitioned_by_weight = partition_binary_numbers_by_weight(n)
        weight = 0
        for i in binary_numbers_partitioned_by_weight:
            if weight != get_hamming_weight(i):
                weight += 1
        assert weight <= length
        sorted_list_of_numbers = sorted(binary_numbers_partitioned_by_weight)
        for i in range(n):
            assert sorted_list_of_numbers[i] == i

    def test_generate_syndrome_decoding_table(self):
        n = 6
        r = 3
        k = n - r
        d = 1
        parity_check_matrix, a_matrix_transposed = fill_parity_check_matrix_and_a_matrix_transposed(n, r, d)
        generator_matrix = fill_generator_matrix_from_a_matrix_transposed(a_matrix_transposed, k, r)
        print()
        print('A parity check matrix:')
        print()
        for num in parity_check_matrix:
            print(format(num, '#0{0}b'.format(n + 2)))
        print()
        print('A generator matrix:')
        print()
        for num in generator_matrix:
            print(format(num, '#0{0}b'.format(n + 2)))
        syndrome_decoding_table = generate_syndrome_decoding_table(generator_matrix, parity_check_matrix, n, k)
        for key in syndrome_decoding_table.keys():
            print()
            print(format(key, '#0{0}b'.format(k + 2)))
            print('---------')
            for word in syndrome_decoding_table[key]:
                print(format(word, '#0{0}b'.format(n + 2)))
        check_list = []
        number_of_words = 2 ** n
        for i in range(number_of_words):
            check_list.append(i)
        number_of_words_in_table = 0
        for key in syndrome_decoding_table.keys():
            for word in syndrome_decoding_table[key]:
                check_list.remove(word)
                number_of_words_in_table += 1
        assert len(check_list) == 0
        assert number_of_words_in_table == number_of_words

    def test_get_binomial_coefficient(self):
        assert get_binomial_coefficient(10, 1) == 10
        assert get_binomial_coefficient(10, 2) == 45
        assert get_binomial_coefficient(10, 3) == 120
        assert get_binomial_coefficient(10, 4) == 210

    def test_calculate_p_error(self):
        n = 10
        t = 1
        p = 0.01
        assert round(calculate_p_error(n, t, p), 5) == 0.00427

    def test_read_file_to_dict(self):
        config = yaml.safe_load(open('config.yml'))
        dictionary = read_file_to_dict(config['decoder-syndrome-decoding'])
        n = 10
        k = 3
        r = n - k
        for key in dictionary.keys():
            print()
            print('{0:0>{width}b}'.format(key, width=r))
            print('---------')
            for word in dictionary[key]:
                print('{0:0>{width}b}'.format(word, width=n))

    def test_code(self):
        config = yaml.safe_load(open('config.yml'))
        for i in range(2 ** self.k):
            code, distorted_code, error = linearcode.code(
                coder_file=config['coder-generator'],
                message=i,
                m_length=self.k,
                error=74)
            assert transpose_matrix(multiply_matrices(
                matrix1=linearcode.read_file_to_list(config['decoder-parity-check']),
                num_col_m1=self.n,
                matrix2=transpose_matrix([code], self.n),
                num_col_m2=1
            ), 1)[0] == 0
            print()
            print('Message:           {0}\n'.format(i))
            print('Code:           {0:0>{width}b}\n'.format(code, width=self.n))
            print('Distorted code: {0:0>{width}b}\n'.format(distorted_code, width=self.n))
            print('Error:          {0:0>{width}b}\n'.format(error, width=self.n))

    def test_decode(self):
        config = yaml.safe_load(open('config.yml'))
        message = linearcode.decode(
            config['decoder-parity-check'],
            self.n,
            config['decoder-syndrome-decoding'],
            int('110010000110111000100', 2)
        )
        print()
        print(message)

    def test_code_decode(self):
        config = yaml.safe_load(open('config.yml'))
        number_of_unfixed_errors = 0
        print()
        number_of_repetitions = 100
        for i in range(number_of_repetitions):
            message = random.randrange(2 ** self.k)
            rand_error = get_random_number_of_hamming_weight(self.n, int((self.d - 1) / 2))
            code, distorted_code, error = linearcode.code(
                coder_file=config['coder-generator'],
                message=message,
                m_length=self.k,
                error=rand_error)
            assert rand_error == error
            decoded = linearcode.decode(
                parity_check_file=config['decoder-parity-check'],
                n=self.n,
                syndrome_file=config['decoder-syndrome-decoding'],
                distorted_code=distorted_code
            )
            if decoded != message:
                number_of_unfixed_errors += 1
                print('Message is {0}, distorted code is {1:0>{width}b}, error is {2:0>{width}b}, deciphered is {3}'
                      .format(message, distorted_code, rand_error, decoded, width=self.n))
        print()
        print('Number of errors which weren\'t fixed is {0}, the rate is {1}'
              .format(number_of_unfixed_errors, number_of_unfixed_errors / number_of_repetitions))

    def test_get_random_number_of_hamming_weight(self):
        for i in range(100):
            length = random.randrange(100)
            if length == 0:
                continue
            weight = random.randrange(length)
            probe = get_random_number_of_hamming_weight(length, weight)
            assert get_hamming_weight(probe) == weight


if __name__ == '__main__':
    unittest.main()
