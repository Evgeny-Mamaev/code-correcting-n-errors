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
        if not is_gilbert_varshamov_bound(n=n, k=k, d=d):
            raise ValueError('The given n == {0}, k == {1} and d == {2} aren\'t compliant'
                             ' with the Gilbert-Varshamov bound.'.format(n, k, d))
        self.n = n
        self.k = k
        self.d = d
        self.r = n - k
        self.parity_check_matrix, self.a_matrix_transposed = fill_parity_check_matrix_and_a_matrix_transposed(n=self.n,
                                                                                                              r=self.r)
        self.generator_matrix = fill_generator_matrix(n=self.n, k=self.k, r=self.r)
        self.i_matrix = fill_i_matrix(size=self.r)


def fill_parity_check_matrix_and_a_matrix_transposed(n, r):
    """
    Produces H-matrix in standard form.
    Fills the 2-d array with the binary
    numbers from 1 to n so, that the
    last columns contain the powers of
    two in decreasing order.

    For example, for n = 7, r = 3:
                0111100
    binary:     1011010
                1101001

    decimal:    3567421

    :param n: number of columns
    :param r: number of rows
    :return: 2-d array filled with binary numbers
    """
    k = n - r
    a_matrix = [0] * r
    powers_of_two, not_powers_of_two = partition_by_is_power_of_two(num=n)
    powers_of_two.reverse()
    not_powers_and_powers_of_two = not_powers_of_two + powers_of_two
    a_matrix_transposed = []
    for i in range(k):
        num = not_powers_and_powers_of_two[i]
        allocate_bits_in_column(
            matrix=a_matrix,
            position=i,
            number=num,
            n=n)
        a_matrix_transposed.append(num)
    return append_matrix_to_front(where=fill_i_matrix(r), what=a_matrix, lshift=0), a_matrix_transposed


def fill_generator_matrix(n, k, r):
    powers_of_two, not_powers_of_two = partition_by_is_power_of_two(num=n)
    not_powers_of_two.reverse()
    powers_and_not_powers_of_two = powers_of_two + not_powers_of_two
    a_matrix = []
    mask_of_r_length = (2 ** r - 1)
    for i in range(k):
        a_matrix.append(powers_and_not_powers_of_two[i] & mask_of_r_length)
    return append_matrix_to_front(where=a_matrix, what=fill_i_matrix(k), lshift=r)


def generator_matrix_from_a_matrix_transposed(transposed_a_matrix, k, r):
    return append_matrix_to_front(
        where=transposed_a_matrix,
        what=fill_i_matrix(k),
        lshift=r)


def transpose_matrix(matrix, number_of_columns):
    transposed_matrix = [0] * number_of_columns
    k = len(matrix)
    for i in range(k):
        num = matrix[i]
        allocate_bits_in_column(
            matrix=transposed_matrix,
            position=i,
            number=num,
            n=k)
    return transposed_matrix


def multiply_matrices(matrix1, num_col_m1, matrix2, num_col_m2):
    m1_length = len(matrix1)
    m2_length = len(matrix2)
    if num_col_m1 != m2_length:
        raise ValueError('The matrices are incompatible: '
                         'matrix1 size is nxm = {0}x{1}, '
                         'matrix2 size is mxp = {2}x{3}. '
                         'm != m '.format(m1_length, num_col_m1, m2_length, num_col_m2))
    transposed_matrix2 = transpose_matrix(matrix2, num_col_m2)
    result = [0] * m1_length
    for i in range(m1_length):
        for j in range(num_col_m2):
            result[i] = result[i] | compute_parity(matrix1[i] & transposed_matrix2[j], m2_length) << (
                        num_col_m2 - j - 1)
    return result


def compute_parity(number, length):
    """
    Computes the parity, which is a number of 1s in
    the binary representation of a given number.
    Can be used to calculate the sum of all 1s.
    :param number: any number to compute the parity.
    :param length: the length of the number,
    specified explicitly.
    :return: 1 if the number of 1s is odd,
    0 if the number of 1s is even.
    """
    number = number & (2 ** length - 1)
    i = 1
    while i <= length:
        number = number ^ (number >> i)
        i *= 2
    return number & 1


def negate_matrix(matrix, r):
    k = len(matrix)
    negated_matrix = [0] * k
    for i in range(k):
        negated_matrix[i] = matrix[i] ^ (2 ** r - 1)
    return negated_matrix


def fill_i_matrix(size):
    i_matrix = []
    for i in range(size):
        i_matrix.append(1 << (size - i - 1))
    return i_matrix


def append_matrix_to_front(where, what, lshift):
    if len(where) != len(what):
        raise ValueError('The given matrices are different in their sizes:'
                         'matrix1.size() == {0}, matrix2.size() == {1}'
                         .format(len(where), len(what)))
    for i in range(len(where)):
        where[i] = where[i] | what[i] << lshift
    return where


def partition_by_is_power_of_two(num):
    powers_of_two = []
    not_powers_of_two = []
    for i in range(num):
        if is_power_of_two(num=i + 1):
            powers_of_two.append(i + 1)
        else:
            not_powers_of_two.append(i + 1)
    return powers_of_two, not_powers_of_two


def allocate_bits_in_column(matrix, position, number, n):
    """

    Not a pure function, which modifies the passed 2-d array by filling
    any given row with the binary representation of the given number

    :param matrix: the matrix is to be filled
    :param position: position of a bit in a binary row, from left to right
    :param number: a number is to be split and allocated in a column
           with the aforementioned position
    :param n: number of columns
    """
    r = len(matrix)
    for j in range(r):
        matrix[j] = matrix[j] | (
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
