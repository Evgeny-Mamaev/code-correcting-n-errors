import random
import unittest
from itertools import combinations

import yaml

import linearcode
from fileutils import remove_files
from linearcode import LinearCode, \
    ErrorCorrectingCode, \
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
    read_file_to_list, \
    read_file_to_dict, \
    sum_modulo_2, \
    transpose_matrix, \
    get_random_number_of_hamming_weight


class LinearCodeTest(unittest.TestCase):
    n = 21
    k = 6
    r = 15
    t = 3
    d = 2 * t + 1

    def test_init_l_c(self):
        LinearCode(n=self.n, k=self.k, d=self.d, channel_error_probability=0.01)

    def test_init_l_c_exception(self):
        with self.assertRaises(ValueError):
            LinearCode(n=3, k=3, d=3, channel_error_probability=0.01)

    def test_init_l_c_timeout(self):
        with self.assertRaises(TimeoutError):
            LinearCode(n=30, k=6, d=6, channel_error_probability=0.01)

    def test_init_e_c_c(self):
        ErrorCorrectingCode(r=self.r, n=self.n, t=self.t, channel_error_probability=0.01)

    def test_is_gilbert_varshamov_bound(self):
        #       /10 - 1\     /10 - 1\             5
        # 1  +  |      |  +  |      |  =  46  >  2  =  32  => false
        #       \   1  /     \   2  /
        assert not is_gilbert_varshamov_bound(n=10, k=5, d=4)
        #       /12 - 1\     /12 - 1\             8
        # 1  +  |      |  +  |      |  =  67  <  2  =  256  => true
        #       \   1  /     \   2  /
        assert is_gilbert_varshamov_bound(n=12, k=4, d=4)

    def test_find_linearly_independent_vector(self):
        assert is_gilbert_varshamov_bound(n=self.n, k=self.k, d=self.d)
        i_matrix = fill_i_matrix(size=self.r)
        for i in range(1000):
            result = find_linearly_independent_vector(
                vector_pool=i_matrix,
                number_of_columns_v_p=self.r,
                linear_combination_size=self.d - 2)
            assert get_hamming_weight(num=result) >= self.d - 1

    def test_f_l_i_v_inner_loop(self):
        linear_combination_size = self.d - 2
        combinations_v_p = []
        vector_pool = fill_i_matrix(size=self.r)
        for i in range(linear_combination_size):
            combinations_v_p.append(list(combinations(vector_pool, i + 1)))
        test_value = 4378  # 1000100011010
        assert get_hamming_weight(num=test_value) < self.d - 1
        for i in range(1000):
            assert f_l_i_v_inner_loop(
                probe=test_value,
                combinations_v_p=combinations_v_p) != test_value

    def test_sum_modulo_2(self):
        test_list = []
        n = 10
        for i in range(n):
            test_list.append(i + 1)
        assert sum_modulo_2(iterable=test_list) == 11

    def test_allocate_bits_in_column(self):
        matrix = [0, 0, 0]
        n = 10
        position = 1
        allocate_bits_in_column(matrix=matrix, position=position, number=2, n=n)
        assert 2 ** (n - 1 - position) == matrix[1]

    def test_fill_parity_check_matrix_and_a_matrix_transposed(self):
        parity_check_matrix, a_matrix_transposed = fill_parity_check_matrix_and_a_matrix_transposed(
            n=self.n,
            r=self.r,
            d=self.d)
        for line in parity_check_matrix:
            assert line != 0
        for line in a_matrix_transposed:
            assert get_hamming_weight(line) >= self.d - 1

    def test_is_power_of_two(self):
        num = 256
        assert is_power_of_two(num=num)
        num = 255
        assert not is_power_of_two(num=num)

    def test_fill_i_matrix(self):
        r = 10
        matrix = fill_i_matrix(size=r)
        i = 0
        for num in matrix:
            assert 2 ** (r - 1 - i) == num
            i += 1

    def test_add_matrix_to_front(self):
        n = 10
        k = 7
        r = n - k
        matrix1 = []
        for i in range(k):
            matrix1.append(i + 1)
        matrix2 = fill_i_matrix(size=k)
        i = 0
        for num in add_matrix_in_front(where=matrix1, what=matrix2, lshift=r):
            # each number consists of the power of two and an element of an I-matrix
            assert num == 2 ** (n - 1 - i) + (i + 1)
            i += 1

    def test_add_matrix_to_front_exception(self):
        with self.assertRaises(ValueError):
            add_matrix_in_front([0], [0, 0], 1)

    def test_generator_matrix_from_a_matrix_transposed(self):
        n = 10
        k = 3
        r = n - k
        d = 1
        parity_check_matrix, a_matrix_transposed = fill_parity_check_matrix_and_a_matrix_transposed(n=n, r=r, d=d)
        for i in fill_generator_matrix_from_a_matrix_transposed(
                a_matrix_transposed=a_matrix_transposed,
                k=k,
                r=r):
            assert get_hamming_weight(i) >= d

    def test_transpose_matrix(self):
        n = 10
        r = 7
        d = 1
        matrix, a_matrix_transposed = fill_parity_check_matrix_and_a_matrix_transposed(n=n, r=r, d=d)
        transposed_matrix = transpose_matrix(matrix=matrix, number_of_columns=n)
        doubly_transposed_matrix = transpose_matrix(matrix=transposed_matrix, number_of_columns=r)
        for i in range(len(doubly_transposed_matrix)):
            assert doubly_transposed_matrix[i] == matrix[i]

    def test_compute_parity(self):
        number_parity_0 = 0b101
        number_parity_exceeding_length_0 = 0b1101
        number_parity_1 = 0b11111111111
        number_parity_exceeding_length_1 = 0b111111111111
        assert compute_parity(number=number_parity_0, length=3) == 0
        assert compute_parity(number=number_parity_exceeding_length_0, length=3) == 0
        assert compute_parity(number=number_parity_1, length=11) == 1
        assert compute_parity(number=number_parity_exceeding_length_1, length=3) == 1

    def test_multiply_matrices(self):
        parity_check_matrix, a_matrix_transposed = fill_parity_check_matrix_and_a_matrix_transposed(
            n=self.n,
            r=self.r,
            d=self.d)
        generator_matrix = fill_generator_matrix_from_a_matrix_transposed(
            a_matrix_transposed=a_matrix_transposed,
            k=self.k,
            r=self.r)
        multiplied_matrix = multiply_matrices(
            matrix1=parity_check_matrix,
            num_col_m1=self.n,
            matrix2=transpose_matrix(generator_matrix, self.n),
            num_col_m2=self.k)
        for i in multiplied_matrix:
            assert i == 0

    def test_multiply_matrices_exception(self):
        with self.assertRaises(ValueError):
            multiply_matrices([0, 0], 2, [0], 1)

    def test_partition_by_is_power_of_two(self):
        n = 100
        powers_of_two, not_powers_of_two = partition_by_is_power_of_two(num=n)
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

    def test_get_hamming_weight(self):
        assert get_hamming_weight(num=7) == 3
        assert get_hamming_weight(num=8) == 1
        assert get_hamming_weight(num=1) == 1

    def test_get_binary_numbers_partitioned_by_weight(self):
        n = 127
        length = len(bin(n)) - 2
        binary_numbers_partitioned_by_weight = partition_binary_numbers_by_weight(num=n)
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
        parity_check_matrix, a_matrix_transposed = fill_parity_check_matrix_and_a_matrix_transposed(n=n, r=r, d=d)
        generator_matrix = fill_generator_matrix_from_a_matrix_transposed(
            a_matrix_transposed=a_matrix_transposed,
            k=k,
            r=r)
        syndrome_decoding_table = generate_syndrome_decoding_table(
            generator_matrix=generator_matrix,
            parity_check_matrix=parity_check_matrix,
            n=n,
            k=k,
            d=d)
        for key in syndrome_decoding_table.keys():
            if key == 0:
                continue
            for word in syndrome_decoding_table[key]:
                assert get_hamming_weight(word) <= self.t
                assert word != 0

    def test_get_binomial_coefficient(self):
        assert get_binomial_coefficient(n=10, k=1) == 10
        assert get_binomial_coefficient(n=10, k=2) == 45
        assert get_binomial_coefficient(n=10, k=3) == 120
        assert get_binomial_coefficient(n=10, k=4) == 210

    def test_calculate_p_error(self):
        n = 10
        t = 1
        p = 0.01
        assert round(calculate_p_error(n=n, t=t, p=p), 5) == 0.00427

    def test_read_file_to_dict(self):
        config = yaml.safe_load(open('config-test.yml'))
        dictionary = read_file_to_dict(name=config['decoder-syndrome-decoding-test'])
        n = 10
        k = 3
        r = n - k
        zero_word_counter = 0
        word_counter = 0
        for key in dictionary.keys():
            for word in dictionary[key]:
                if word == 0:
                    zero_word_counter += 1
                word_counter += 1
        assert zero_word_counter <= 1
        assert word_counter > 2 ** r

    def test_code(self):
        config = yaml.safe_load(open('config-test.yml'))
        for i in range(2 ** self.k):
            code, distorted_code, error = linearcode.encode(
                coder_file=config['coder-generator-test'],
                message=i,
                m_length=self.k,
                error=74)  # 000000000000001001010
            assert transpose_matrix(
                matrix=multiply_matrices(
                    matrix1=read_file_to_list(config['decoder-parity-check-test']),
                    num_col_m1=self.n,
                    matrix2=transpose_matrix([code], self.n),
                    num_col_m2=1
                ),
                number_of_columns=1)[0] == 0
        for i in range(2 ** self.k):
            code, distorted_code, error = linearcode.encode(
                coder_file=config['coder-generator-test'],
                message=i,
                m_length=self.k)
            if get_hamming_weight(error) <= self.t:
                assert transpose_matrix(
                    matrix=multiply_matrices(
                        matrix1=read_file_to_list(config['decoder-parity-check-test']),
                        num_col_m1=self.n,
                        matrix2=transpose_matrix([code], self.n),
                        num_col_m2=1
                    ),
                    number_of_columns=1)[0] == 0

    def test_code_exception(self):
        config = yaml.safe_load(open('config-test.yml'))
        with self.assertRaises(ValueError):
            linearcode.encode(
                coder_file=config['coder-generator-test'],
                message=1,
                m_length=2,
                error=74)

    def test_decode(self):
        message = 56
        config = yaml.safe_load(open('config-test.yml'))
        decoded_message = linearcode.decode(
            parity_check_file=config['decoder-parity-check-test'],
            n=self.n,
            syndrome_file=config['decoder-syndrome-decoding-test'],
            distorted_code=int('110010000110111000100', 2)
        )
        assert decoded_message == message

    def test_code_decode(self):
        config = yaml.safe_load(open('config-test.yml'))
        number_of_unfixed_errors = 0
        number_of_repetitions = 100
        for i in range(number_of_repetitions):
            message = random.randrange(2 ** self.k)
            rand_error = get_random_number_of_hamming_weight(self.n, int((self.d - 1) / 2))
            code, distorted_code, error = linearcode.encode(
                coder_file=config['coder-generator-test'],
                message=message,
                m_length=self.k,
                error=rand_error)
            assert rand_error == error
            decoded = linearcode.decode(
                parity_check_file=config['decoder-parity-check-test'],
                n=self.n,
                syndrome_file=config['decoder-syndrome-decoding-test'],
                distorted_code=distorted_code
            )
            if decoded != message:
                number_of_unfixed_errors += 1
        assert number_of_unfixed_errors == 0

    def test_get_random_number_of_hamming_weight(self):
        for i in range(100):
            length = random.randrange(100)
            if length == 0:
                continue
            weight = random.randrange(length)
            probe = get_random_number_of_hamming_weight(length=length, weight=weight)
            assert get_hamming_weight(num=probe) == weight

    def test_get_random_number_of_hamming_weight_exception(self):
        with self.assertRaises(ValueError):
            get_random_number_of_hamming_weight(length=2, weight=3)

    @classmethod
    def tearDownClass(cls) -> None:
        remove_files()


if __name__ == '__main__':
    unittest.main()
