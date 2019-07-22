import array

import math


class LinearCode(object):
    """
    Attributes:
        i_matrix: An (n - k) * (n - k) unit matrix, which looks like
            [ 1 ... 0]
            [ 0 1 ...]
            [ 0 0 1  ]
            [ 0 ... 1].
        a_matrix:    d.
    """

    def __init__(self, n, k, d):
        """

        :param n: codeword length
        :param k: number of significant bits
        :param d: minimum code distance
        """
        self.n = n
        self.k = k
        self.d = d
        self.r = n - k
        self.generator_matrix = array.array('L')  # use 32 bit unsigned integer
        self.parity_check_matrix = [self.r]
        self.a_matrix = []

        self.i_matrix = []
        # fill the unit matrix
        for i in range(n - k):
            self.i_matrix.append(1 << (n - k - i - 1))

        if not is_gilbert_varshamov_bound(n=n, k=k, d=d):
            raise ValueError('The given n == {0}, k == {1} and d == {2} aren\'t compliant'
                             ' with the Gilbert-Varshamov bound.'.format(n, k, d))


def fill_parity_check_matrix(n, r):
    """
    Fills the 2-d array with the binary numbers from 1 to n,
    column 0 contains binary 1 ...
    column n - 1 contains binary n.

    For example, for n = 7, r = 3:
                0001111
    binary:     0110011
                1010101

    decimal:    1234567

    :param n: number of columns
    :param r: number of rows
    :return: 2-d array filled with binary numbers
    """
    parity_check_matrix = [0] * r
    powers_of_two = []
    not_powers_of_two = []
    for i in range(n):
        if is_power_of_two(num=i + 1):
            powers_of_two.append(i + 1)
        else:
            not_powers_of_two.append(i + 1)
    i = 0
    for num in not_powers_of_two:
        allocate_bits_in_column(parity_check_matrix=parity_check_matrix, position=i, number=num, n=n)
        i += 1
    powers_of_two.reverse()
    for num in powers_of_two:
        allocate_bits_in_column(parity_check_matrix=parity_check_matrix, position=i, number=num, n=n)
        i += 1
    return parity_check_matrix


def allocate_bits_in_column(parity_check_matrix, position, number, n):
    """

    Not a pure function, which modifies the passed 2-d array by filling
    any given row with the binary representation of the given number

    :param parity_check_matrix: the matrix is to be filled
    :param position: position of a bit in a binary row, from left to right
    :param number: a number is to be split and allocated in a column
           with the aforementioned position
    :param n: number of columns
    """
    r = len(parity_check_matrix)
    for j in range(r):
        parity_check_matrix[j] = parity_check_matrix[j] | (
            (1 << (n - 1 - position)) if number & (1 << (r - 1 - j)) else 0)


def is_power_of_two(num):
    """
    Determines the binary length of the number and
    checks whether it has only one 1 => the power of 2
    :param num: a number to check
    :return: true if the number has only one 1
             and all the remaining 0,
             false otherwise
    """
    counter = 0
    for i in (range(len(bin(num)) - 2)):
        if num & 1:
            counter += 1
        num >>= 1
        i += 1
    return counter == 1


def is_gilbert_varshamov_bound(n, k, d):
    """

    :param n: codeword length
    :param k: number of significant bits
    :param d: minimum code distance
    :return: the result of the Gilbert-Varshamov bound,
             true if the tuple (n, k, d) satisfies the bound,
             false if not

                /n - 1\              /n - 1\      r
         1  +  |      |  +  ...  +  |      |  <  2 , r = n - k
               \  1  /              \d - 2/
    """
    result = 1
    for i in range(1, d - 1):
        result += math.factorial(n - 1) / math.factorial(i) / math.factorial(n - 1 - i)
    is_in_bound = result < 2 ** (n - k)
    return is_in_bound


class ErrorCorrectingCode(LinearCode):
    def __init__(self, r, n, t):
        """

        :param r: number of check bits
        :param n: codeword length
        :param t: number of errors, which the code should correct
        """
        LinearCode.__init__(self, n, n - r, 2 * t + 1)
