import deprecation
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

        :param n: codeword length.
        :param k: number of significant bits.
        :param d: minimum code distance.
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
        self.generator_matrix = fill_generator_matrix_from_a_matrix_transposed(self.a_matrix_transposed, self.k, self.r)
        self.i_matrix = fill_i_matrix(size=self.r)


def fill_parity_check_matrix_and_a_matrix_transposed(n, r):
    """
    Produces an H-matrix in standard form
    and a transposed A-matrix, corresponding
    to the H-matrix. Fills the H-matrix with
    the binary numbers from 1 to n so, that
    the last columns contain the powers of
    two in decreasing order, and the A-matrix
    contains only non-powers of two.

    For example, for n = 7, r = 3:
                0111100
    binary:     1011010
                1101001
                ____
                 A
                _______
                   H

    decimal:    3567421

    :param n: a number of columns.
    :param r: a number of rows.
    :return: a tuple of two 2-d arrays
    filled with binary numbers:
    the H-matrix in standard form and
    the A-matrix transposed.
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
    return add_matrix_in_front(where=fill_i_matrix(r), what=a_matrix, lshift=0), a_matrix_transposed


def fill_generator_matrix_from_a_matrix_transposed(a_matrix_transposed, k, r):
    """
    Fills a G-matrix from a transposed A-matrix.
        [    |  tr]
    G = [ I  | A  ]
        [  k |    ]
                                    tr
    :param a_matrix_transposed: an A  -matrix.
    :param k: the number of significant bits.
    :param r: the number of check bits.
    :return: the G-matrix.
    """
    return add_matrix_in_front(
        where=a_matrix_transposed,
        what=fill_i_matrix(k),
        lshift=r)


def transpose_matrix(matrix, number_of_columns):
    """
    Transposes a given matrix.
    :param matrix: a matrix to transpose.
    :param number_of_columns: the length of
    the numbers in the given matrix, the same
    as the number of columns in the matrix.
    :return: a transposed matrix.
    """
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
    """
    Multiplies two matrices with the given number of columns.
    :param matrix1: the first matrix to multiply.
    :param num_col_m1: the number of columns in
    the first matrix.
    :param matrix2: the second matrix to multiply.
    :param num_col_m2: the number of columns in
    the second matrix.
    :raises: ValueError in case if the number of
    columns of the first matrix doesn't equal to
    the number of lines of the second matrix.
    :return: the product of two matrices.
    """
    m1_length = len(matrix1)
    m2_length = len(matrix2)
    if num_col_m1 != m2_length:
        raise ValueError('The matrices are incompatible: '
                         'matrix1 size is nxm = {0}x{1}, '
                         'matrix2 size is mxp = {2}x{3}. '
                         '{1} != {2} '.format(m1_length, num_col_m1, m2_length, num_col_m2))
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


@deprecation.deprecated(deprecated_in="1.0", removed_in="2.0",
                        current_version="1.0",
                        details="An unused function")
def negate_matrix(matrix, r):
    k = len(matrix)
    negated_matrix = [0] * k
    for i in range(k):
        negated_matrix[i] = matrix[i] ^ (2 ** r - 1)
    return negated_matrix


def fill_i_matrix(size):
    """
    Fills an I-matrix (unit):
    [1 ... 0]
    [ 1... 0]
    [  ...  ]
    [0 ... 1]
    :param size: the size of the unit matrix.
    :return: a unit matrix.
    """
    i_matrix = []
    for i in range(size):
        i_matrix.append(1 << (size - i - 1))
    return i_matrix


def add_matrix_in_front(where, what, lshift):
    """
    Adds one matrix in front of another with
    the provided shift, eg:

    [100]                                   [11]
    [010] is to be appended to the front of [00] =>
    [001]                                   [11]
    what                                    where
        [10011]
     => [01000] , lshift = 2.
        [00111]

    :param where: the matrix where to add another one.
    :param what: the matrix which to add.
    :param lshift: the number of columns to shift,
    normally it's the width of the where matrix.
    :raises: ValueError if the numbers of lines
    of the matrices are different.
    :return: the sum of two matrices in terms of
    extending each line of one matrix by the lines
    of another matrix.
    """
    if len(where) != len(what):
        raise ValueError('The given matrices are different in their'
                         ' number of lines:'
                         ' len(matrix1) == {0}, len(matrix2) == {1}'
                         .format(len(where), len(what)))
    for i in range(len(where)):
        where[i] = where[i] | what[i] << lshift
    return where


def partition_by_is_power_of_two(num):
    """
    Partitions the binary representation
    of all the preceding numbers including
    the given number by the criteria is it
    the power of 2.
    :param num: a number to analyse.
    :return: a tuple of two lists:
    the first one is the powers of 2,
    the second one is not the powers of 2.
    """
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
    checks whether it has only one 1 => the power of 2.
    :param num: a number to check.
    :return: true if the number has only one 1
             and all the remaining 0,
             false otherwise.
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
    Computes the Gilbert-Varshamov bound and
    determines whether the specified combination
    of parameters complies with the bound.
    :param n: codeword length
    :param k: number of significant bits
    :param d: minimum code distance
    :return: the result of the Gilbert-Varshamov bound,
             true if the tuple (n, k, d) satisfies the bound,
             false if not.

                /n - 1\              /n - 1\      r
         1  +  |      |  +  ...  +  |      |  <  2 , r = n - k
               \  1  /              \d - 2/
    """
    result = 1
    for i in range(1, d - 1):
        result += math.factorial(n - 1) / math.factorial(i) / math.factorial(n - 1 - i)
    is_in_bound = result < 2 ** (n - k)
    return is_in_bound


def generate_syndrome_decoding_table(generator_matrix, parity_check_matrix, n, k):
    """
    Generates a syndrome decoding table,
    which maps a coset to its syndrome S.
    Similar to a standard array, contains
    the same data, but is arranged differently.
    Let L be an [n, k] linear code. For any
    vector a, the set a + L = {a + x:x E L}
    is called a coset (or translate) of L.
             tr
    S = H * y  , e.g. mapping is as follows:

    S:    message:   00    10    01    11
    [0]
    [0]  ->  code:  0000  1011  0101  1110

    [1]
    [1]  -> coset:  1000  0011  1101  0110

    [0]
    [1]  -> coset:  0100  1111  0001  1010

    [1]
    [0]  -> coset:  0010  1001  0111  1100
                   coset
                  leaders
    :param generator_matrix: a generator matrix
    :param parity_check_matrix: a parity check
    matrix.
    :param n: the length of the code word
    :param k: the number of significant bits
    :return: a table, which maps lists of error
    vectors belong to the same coset to the
    coset syndrome.
    """
    messages = 2 ** k
    r = n - k
    cosets = 2 ** r
    coset_shifts = 2 ** n
    syndrome_decoding_table = {}
    for i in range(2 ** r):
        syndrome_decoding_table[i] = []
    binary_numbers_partitioned_by_weight = partition_binary_numbers_by_weight(coset_shifts)
    counter = 0
    for j in range(coset_shifts):
        if counter == cosets:
            break
        coset_shift = binary_numbers_partitioned_by_weight[j]
        coset_list = []
        for i in range(messages):
            code = multiply_matrices([i], k, generator_matrix, n)[0]
            coset_vector = code ^ coset_shift
            if get_hamming_weight(coset_vector) < get_hamming_weight(coset_shift):
                coset_list = []
                break
            coset_list.append(coset_vector)
        if len(coset_list) == 0:
            continue
        syndrome_transposed = multiply_matrices(parity_check_matrix, n, transpose_matrix([counter], n), 1)
        syndrome = transpose_matrix(syndrome_transposed, 1)[0]
        for i in coset_list:
            syndrome_decoding_table[syndrome].append(i)
        counter += 1
    return syndrome_decoding_table


def partition_binary_numbers_by_weight(num):
    """
    Partitions all the numbers until
    the specified number by the weight
    of their binary representations.
    :param num: a number until which
    all the numbers are to be
    partitioned.
    :return: a list of partitioned
    by the weight numbers, increasing
    order.
    """
    length = len(bin(num)) - 2
    weight_map = {}
    for i in range(length):
        weight_map[i] = []
    for i in range(num):
        weight = get_hamming_weight(i)
        weight_map[weight].append(i)
    result = []
    for list_of_numbers in weight_map.values():
        for i in list_of_numbers:
            result.append(i)
    return result


def get_hamming_weight(num):
    """
    Calculates the Hamming weight.
    [10001] the weight is 2,
    [10] the weight is 1.
    :param num: a number to calculate
    the weight of.
    :return: the Hamming weight of
    the binary representation of
    the given number.
    """
    length = len(bin(num)) - 2
    weight = 0
    for i in range(length):
        if num & 1 == 1:
            weight += 1
        num >>= 1
    return weight


class ErrorCorrectingCode(LinearCode):
    def __init__(self, r, n, t):
        """

        :param r: number of check bits
        :param n: codeword length
        :param t: number of errors, which the code should correct
        """
        LinearCode.__init__(self, n, n - r, 2 * t + 1)
