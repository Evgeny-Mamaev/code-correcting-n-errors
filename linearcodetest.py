import unittest

from linearcode import LinearCode
from linearcode import allocate_bits_in_column
from linearcode import append_matrix_to_front
from linearcode import compute_parity
from linearcode import fill_generator_matrix
from linearcode import fill_i_matrix
from linearcode import fill_parity_check_matrix_and_a_matrix_transposed
from linearcode import generator_matrix_from_a_matrix_transposed
from linearcode import is_gilbert_varshamov_bound
from linearcode import is_power_of_two
from linearcode import multiply_matrices
from linearcode import negate_matrix
from linearcode import partition_by_is_power_of_two
from linearcode import transpose_matrix

# TODO test exceptions
class LinearCodeTest(unittest.TestCase):
    def test_init(self):
        linear_code = LinearCode(n=10, k=3, d=3)
        print(bin(linear_code.i_matrix[1]))
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

    def test_allocate_bits_in_column(self):
        matrix = [0, 0, 0]
        n = 10
        position = 1
        allocate_bits_in_column(matrix=matrix, position=position, number=2, n=n)
        assert 2 ** (n - 1 - position) == matrix[1]
        pass

    def test_fill_parity_check_matrix_and_a_matrix_transposed(self):
        n = 10
        r = 7
        parity_check_matrix, a_matrix_transposed = fill_parity_check_matrix_and_a_matrix_transposed(n, r)
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
        for num in append_matrix_to_front(matrix1, matrix2, r):
            # each number consists of the power of two and an element of an I-matrix
            assert num == 2 ** (n - 1 - i) + (i + 1)
            i += 1
        pass

    def test_fill_generator_matrix(self):
        n = 10
        k = 7
        r = n - k
        print()
        for i in fill_generator_matrix(n, k, r):
            print(format(i, '#0{0}b'.format(n + 2)))
        pass

    def test_negate_matrix(self):
        n = 5
        print()
        for i in negate_matrix(fill_i_matrix(n), 10):
            print(format(i, '#0{0}b'.format(n + 2)))
        pass

    def test_generator_matrix_from_a_matrix_transposed(self):
        n = 10
        k = 3
        r = n - k
        parity_check_matrix, a_matrix_transposed = fill_parity_check_matrix_and_a_matrix_transposed(n, r)
        print()
        for i in parity_check_matrix:
            print(format(i, '#0{0}b'.format(n + 2)))
        print()
        for i in generator_matrix_from_a_matrix_transposed(a_matrix_transposed, k, r):
            print(format(i, '#0{0}b'.format(n + 2)))

    def test_transpose_matrix(self):
        n = 10
        r = 7
        matrix, a_matrix_transposed = fill_parity_check_matrix_and_a_matrix_transposed(n, r)
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
        n = 10
        k = 3
        r = n - k
        parity_check_matrix, a_matrix_transposed = fill_parity_check_matrix_and_a_matrix_transposed(n, r)
        generator_matrix = generator_matrix_from_a_matrix_transposed(a_matrix_transposed, k, r)
        multiplied_matrix = multiply_matrices(parity_check_matrix, n, transpose_matrix(generator_matrix, n), k)
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


if __name__ == '__main__':
    unittest.main()
