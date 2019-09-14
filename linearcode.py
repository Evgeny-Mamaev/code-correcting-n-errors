import random

import math
import yaml
from itertools import combinations

import timeout


class LinearCode(object):

    @timeout.timeout(yaml.safe_load(open("config.yml"))['timeout'])
    def __init__(self, n, k, d, channel_error_probability=0.01):
        """

        :param n: a codeword length.
        :param k: a number of significant bits.
        :param d: a minimum code distance.
        :param channel_error_probability: a
        channel error probability.
        :raises ValueError if the specified
        parameters aren't complaint with the
        Gilbert-Varshamov bound.
        """
        if not is_gilbert_varshamov_bound(n=n, k=k, d=d):
            raise ValueError('The given n == {0}, k == {1} and d == {2} aren\'t compliant'
                             ' with the Gilbert-Varshamov bound.'.format(n, k, d))
        self.channel_error_probability = channel_error_probability
        self.n = n
        self.k = k
        self.d = d
        self.r = n - k
        self.parity_check_matrix, self.a_matrix_transposed = fill_parity_check_matrix_and_a_matrix_transposed(
            n=self.n,
            r=self.r,
            d=self.d)
        self.generator_matrix = fill_generator_matrix_from_a_matrix_transposed(
            a_matrix_transposed=self.a_matrix_transposed,
            k=self.k,
            r=self.r)
        self.syndrome_decoding_table = generate_syndrome_decoding_table(
            generator_matrix=self.generator_matrix,
            parity_check_matrix=self.parity_check_matrix,
            n=n,
            k=k,
            d=d)
        self.i_matrix = fill_i_matrix(size=self.r)
        print_files(
            generator_matrix=self.generator_matrix,
            parity_check_matrix=self.parity_check_matrix,
            syndrome_decoding_table=self.syndrome_decoding_table,
            n=n,
            k=k)
        print_code_error_probability(
            n=n,
            t=((d - 1) / 2).as_integer_ratio()[0],
            p=channel_error_probability)


def fill_parity_check_matrix_and_a_matrix_transposed(n, r, d):
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
    :param d: a code distance.
    :return: a tuple of two 2-d arrays
    filled with binary numbers:
    the H-matrix in standard form and
    the A-matrix transposed.
    """
    k = n - r
    a_matrix = [0] * r
    vector_pool = fill_i_matrix(r)
    a_matrix_transposed = []
    for i in range(k):
        num = find_linearly_independent_vector(
            vector_pool=vector_pool,
            number_of_columns_v_p=r,
            linear_combination_size=d - 2)
        allocate_bits_in_column(
            matrix=a_matrix,
            position=i,
            number=num,
            n=n)
        vector_pool.insert(i, num)
        a_matrix_transposed.append(num)
    return add_matrix_in_front(where=fill_i_matrix(r), what=a_matrix, lshift=0), a_matrix_transposed


def find_linearly_independent_vector(vector_pool, number_of_columns_v_p, linear_combination_size):
    """
    Finds a random number which meets the criteria
    of linear independence in terms of the binary
    representation from all the possible linear 
    combinations of size less then or equal to the 
    specified size of the specified pool of numbers.

    How the algorithms works:
    It is a stochastic process of finding a linearly
    independent vector. The algorithm iterates through
         / i \    / i \      / i \
    all |    | , |    | ... |    | linear combinations
        \ 1 /    \ 2 /      \ n /
    of vectors from
    :param vector_pool, where i stands for the number
    of vectors in the pool and n stands for
    :param linear_combination_size to check whether
    the randomly picked vector is linearly independent
    from the combinations. For further information refer
    to the Gilbert-Varshamov bound theorem proof.
    :param vector_pool: a list of numbers which
    the function finds a linearly independent number
    from in terms of its binary representation.
    :param number_of_columns_v_p: a length of
    numbers in the pool.
    :param linear_combination_size: a size of
    the set which forms a liner combination.
    :return: a linearly independent binary vector
    in form of int.
    """
    combinations_v_p = []
    for i in range(linear_combination_size):
        combinations_v_p.append(list(combinations(vector_pool, i + 1)))
    while True:
        probe = random.randrange(2 ** number_of_columns_v_p)
        result = f_l_i_v_inner_loop(probe=probe, combinations_v_p=combinations_v_p)
        if result != 0:
            return probe


def f_l_i_v_inner_loop(probe, combinations_v_p):
    """
    Participates in the search of the linearly
    independent vector as an inner loop.
    :param probe: a random vector to determine
    its linear independence.
    :param combinations_v_p: combinations
    of vectors, which form linear combinations.
    The sum modulo 2 of each of these
    combinations is examined to determine the
    linear dependence of the probe. If none of
    them match the probe, it's linearly
    independent.
    :return: 0 if the probe is equal to any of
    the sum, the probe itself otherwise.
    """
    for particular_combinations in combinations_v_p:
        for combination in particular_combinations:
            if probe == sum_modulo_2(iterable=combination):
                return 0
    return probe


def sum_modulo_2(iterable):
    """
    Sums modulo 2 elements of an iterable.
    :param iterable: an iterable.
    :return: the sum module 2 of the elements.
    """
    result = 0
    for i in iterable:
        result = result ^ i
    return result


def get_random_number_of_hamming_weight(length, weight):
    """
    Returns a random number of a particular Hamming
    weight.
    :param length: a length of a number.
    :param weight: a weight of a number.
    :return: the number meets the requirements of
    length and Hamming weight.
    :raises: ValueError if the weight is greater
    than the length.
    """
    if weight > length:
        raise ValueError('The weight shouldn\'t be greater'
                         ' than the length: {0} > {1}'
                         .format(weight, length))
    i = 0
    result = 0
    while True:
        if i == weight:
            return result
        shift = random.randrange(length)
        power_of_two = 1 << shift
        if power_of_two & result == power_of_two:
            continue
        result |= power_of_two
        i += 1


def fill_generator_matrix_from_a_matrix_transposed(a_matrix_transposed, k, r):
    """
    Fills a G-matrix from a transposed A-matrix.
        [    |  tr]
    G = [ I  | A  ]
        [  k |    ]
                                    tr
    :param a_matrix_transposed: an A  -matrix.
    :param k: a number of significant bits.
    :param r: a number of check bits.
    :return: a G-matrix.
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
    :return: the product of two matrices.
    :raises: ValueError in case if the number of
    columns of the first matrix doesn't equal to
    the number of lines of the second matrix.
    """
    m1_length = len(matrix1)
    m2_length = len(matrix2)
    if num_col_m1 != m2_length:
        raise ValueError('The matrices are incompatible: '
                         'matrix1 size is nxm = {0}x{1}, '
                         'matrix2 size is mxp = {2}x{3}. '
                         '{1} != {2} '.format(m1_length, num_col_m1, m2_length, num_col_m2))
    transposed_matrix2 = transpose_matrix(matrix=matrix2, number_of_columns=num_col_m2)
    result = [0] * m1_length
    for i in range(m1_length):
        for j in range(num_col_m2):
            result[i] = result[i] | compute_parity(
                number=matrix1[i] & transposed_matrix2[j],
                length=m2_length) << (num_col_m2 - j - 1)
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
    :param lshift: a number of columns to shift,
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
    any given row with the binary representation of the given number.
    :param matrix: the matrix is to be filled.
    :param position: a position of a bit in a binary row, from left to right.
    :param number: a number is to be split and allocated in a column
    with the aforementioned position.
    :param n: number of columns.
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
    and all the remaining 0, false otherwise.
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
    :param n: a codeword length
    :param k: a number of significant bits
    :param d: a minimum code distance
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


def generate_syndrome_decoding_table(generator_matrix, parity_check_matrix, n, k, d):
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
    [1]  -> coset:  1000

    [0]
    [1]  -> coset:  0100

    [1]
    [0]  -> coset:  0010
                   coset
                  leaders
    :param generator_matrix: a generator matrix.
    :param parity_check_matrix: a parity check
    matrix.
    :param n: the length of the code word.
    :param k: a number of significant bits.
    :param d: a code distance.
    :return: a table, which maps lists of error.
    vectors have the same syndrome to the syndrome.
    """
    messages = 2 ** k
    r = n - k
    cosets = 2 ** r
    coset_shifts = 2 ** n
    syndrome_decoding_table = {0: []}
    for i in range(messages):
        syndrome_decoding_table[0].insert(i, multiply_matrices(
            matrix1=[i],
            num_col_m1=k,
            matrix2=generator_matrix,
            num_col_m2=n)[0])
    binary_numbers_partitioned_by_weight = partition_binary_numbers_by_weight(num=coset_shifts)
    j = 0
    for i in binary_numbers_partitioned_by_weight:
        syndrome_transposed = multiply_matrices(
            matrix1=parity_check_matrix,
            num_col_m1=n,
            matrix2=transpose_matrix(
                matrix=[i],
                number_of_columns=n),
            num_col_m2=1)
        syndrome = transpose_matrix(
            matrix=syndrome_transposed,
            number_of_columns=1)[0]
        if syndrome == 0:
            continue
        if get_hamming_weight(i) <= (d - 1) / 2:
            if syndrome in syndrome_decoding_table:
                syndrome_decoding_table[syndrome].append(i)
            else:
                syndrome_decoding_table[syndrome] = []
                syndrome_decoding_table[syndrome].append(i)
        if j >= cosets:
            break
        j += 1
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
        weight = get_hamming_weight(num=i)
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


def get_binomial_coefficient(n, k):
    """
    Calculates a binomial coefficient
     / n \
    |    |
    \ k /.
    :param n: a number of elements
    from which the combinations should
    be picked.
    :param k: a number of elements in
    each combination.
    :return: a number of combinations.
    """
    r = n - k
    return math.factorial(n) / math.factorial(k) / math.factorial(r)


def calculate_p_error(n, t, p):
    """
    Calculates the probability of
    error for the code with the
    particular parameters.
    :param n: a codeword length.
    :param t: a number of correctable
    errors.
    :param p: the probability of
    the channel error.
    :return: the probability of
    error for the code correcting
    the specified number of errors.
    """
    sum = 0
    for i in range(t + 1):
        sum += get_binomial_coefficient(n=n, k=i) * p ** i * (1 - p) ** (n - i)
    return 1 - sum


def print_files(generator_matrix, parity_check_matrix, syndrome_decoding_table, n, k):
    """
    Prints files for the coding-
    decoding process.
    :param generator_matrix: a generator
    matrix used for a coder.
    :param parity_check_matrix: a parity
    check matrix used for a decoder.
    :param syndrome_decoding_table:
    a syndrome decoding table used for
    a decoder to match a distorted vector
    with a coset leader and then find
    the difference of the distorted vector
    and the leader among the vectors from
    the zeroth coset.
    :param n: a length of a codeword.
    :param k: a length of a message.
    """
    config = yaml.safe_load(open("config.yml"))
    r = n - k
    with open(config['coder-generator'], 'a') as file:
        for i in generator_matrix:
            file.write('{0:0>{width}b}'.format(i, width=n))
            file.write('\n')
    with open(config['decoder-parity-check'], 'a') as file:
        for i in parity_check_matrix:
            file.write('{0:0>{width}b}'.format(i, width=n))
            file.write('\n')
    with open(config['decoder-syndrome-decoding'], 'a') as file:
        for key in syndrome_decoding_table.keys():
            file.write('\n')
            file.write('S:{0:0>{width}b}'.format(key, width=r))
            file.write('\n')
            file.write('---------')
            file.write('\n')
            for word in syndrome_decoding_table[key]:
                file.write('{0:0>{width}b}'.format(word, width=n))
                file.write('\n')


def print_code_error_probability(n, t, p):
    """
    Prints the error probability of a code.
    :param n: a codeword length.
    :param t: a number of correctable
    errors.
    :param p: the probability of
    the channel error.
    """
    print()
    print('The decoder error probability is {0:.2%}'.format(calculate_p_error(n, t, p)))


def read_file_to_list(name):
    """
    Reads a file consists of a
    list of binary vectors.
    :param name: a name of
    a file.
    :return: a list of binary
    vectors in form of int.
    """
    lines = []
    with open(name) as file:
        for line in file:
            line = line.strip()
            lines.append(int(line, 2))
    return lines


def read_file_to_dict(name):
    """
    Reads a file consists of a
    dictionary in a special form:
    S:000
    ___________
    0000001
    0000010
    where S: stands for a syndrome
    and the vectors below the dash
    stand for the error vectors,
    which have the particular
    syndrome. A syndrome is a key,
    the vectors are the values of
    a dictionary.
    :param name: a name of
    a file.
    :return: a dict of binary
    vectors map to a list of
    error vectors in form of int.
    """
    dictionary = {}
    with open(name) as file:
        key = 0
        for line in file:
            line = line.rstrip()
            if line.isdigit():
                dictionary[key].append(int(line, 2))
            if line.startswith('S'):
                key = int(line.partition(':')[2], 2)
                dictionary[key] = []
    return dictionary


@timeout.timeout(yaml.safe_load(open("config.yml"))['timeout'])
def encode(coder_file, message, m_length, error=0):
    """
    Codes a message and distorts the code
    by the specified error vector if any
    presents or an arbitrary error.
    :param coder_file: a file where a
    generator matrix is written.
    :param message: a message to code.
    :param m_length: a message length.
    :param error: an optional error,
    which distorts an obtained code.
    :return: a tuple of a code, a
    distorted code, an error vector.
    :raises: ValueError if the length
    of the G-matrix and m_length are
    incompatible.
    """
    generator_matrix = read_file_to_list(name=coder_file)
    max_line = max(generator_matrix)
    number_of_columns_g_m = len(bin(max_line)) - 2
    gen_matrix_length = len(generator_matrix)
    if gen_matrix_length != m_length:
        raise ValueError('The given len(generator_matrix) == {0}'
                         ' and m_length == {1} are incompatible.'
                         ' {0} != {1}.'.
                         format(gen_matrix_length, m_length))
    if error == 0:
        error = random.randrange(2 ** number_of_columns_g_m)
    code = multiply_matrices(
        matrix1=[message],
        num_col_m1=m_length,
        matrix2=generator_matrix,
        num_col_m2=number_of_columns_g_m)[0]
    return code, code ^ error, error


@timeout.timeout(yaml.safe_load(open("config.yml"))['timeout'])
def decode(parity_check_file, n, syndrome_file, distorted_code):
    """
    Decodes a distorted message.
    :param parity_check_file: a file
    contains a parity check matrix.
    :param n: a number of columns in
    the parity check matrix.
    :param syndrome_file: a file
    contains a syndrome decoding
    table.
    :param distorted_code: a
    distorted code.
    :return: a decoded message
    """
    parity_check_matrix = read_file_to_list(name=parity_check_file)
    syndrome_decoding_table = read_file_to_dict(name=syndrome_file)
    syndrome = transpose_matrix(
        matrix=multiply_matrices(
            matrix1=parity_check_matrix,
            num_col_m1=n,
            matrix2=transpose_matrix(
                matrix=[distorted_code],
                number_of_columns=n),
            num_col_m2=1),
        number_of_columns=1)[0]
    return syndrome_decoding_table[0].index(distorted_code ^ syndrome_decoding_table[syndrome][0])


class ErrorCorrectingCode(LinearCode):
    @timeout.timeout(yaml.safe_load(open("config.yml"))['timeout'])
    def __init__(self, r, n, t, channel_error_probability=0.01):
        """

        :param r: a number of check bits.
        :param n: a codeword length.
        :param t: a number of errors, which the code should correct.
        :param channel_error_probability: a
        channel error probability.
        """
        LinearCode.__init__(
            self=self,
            n=n,
            k=n - r,
            d=2 * t + 1,
            channel_error_probability=channel_error_probability)
