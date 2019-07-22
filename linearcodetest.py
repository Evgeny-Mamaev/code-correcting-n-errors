import unittest

from linearcode import LinearCode
from linearcode import allocate_bits_in_column
from linearcode import fill_parity_check_matrix
from linearcode import is_gilbert_varshamov_bound
from linearcode import is_power_of_two


class LinearCodeTest(unittest.TestCase):
    def test_init(self):
        linear_code = LinearCode(10, 5, 4)
        print(bin(linear_code.i_matrix[1]))
        assert True
        pass

    def test_gilbert_varshamov_false(self):
        #       /10 - 1\     /10 - 1\             5
        # 1  +  |      |  +  |      |  =  46  >  2  =  32  => false
        #       \   1  /     \   2  /
        assert not is_gilbert_varshamov_bound(10, 5, 4)
        pass

    def test_gilbert_varshamov_true(self):
        #       /12 - 1\     /12 - 1\             8
        # 1  +  |      |  +  |      |  =  67  <  2  =  256  => true
        #       \   1  /     \   2  /
        assert is_gilbert_varshamov_bound(12, 4, 4)
        pass

    def test_allocate_bits_in_column(self):
        matrix = [0, 0, 0]
        n = 10
        position = 1
        allocate_bits_in_column(parity_check_matrix=matrix, position=position, number=2, n=n)
        assert 2 ** (n - 1 - position) == matrix[1]
        pass

    def test_fill_parity_check_matrix(self):
        n = 10
        r = 3
        parity_check_matrix = fill_parity_check_matrix(n, r)
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


if __name__ == '__main__':
    unittest.main()
