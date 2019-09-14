import yaml

from fileutils import remove_files
from linearcode import ErrorCorrectingCode, encode, decode, get_random_number_of_hamming_weight


def simulate():
    remove_files()
    r = int(input("Enter number of check bits r: "))
    n = int(input("Enter length of a code word n: "))
    t = int(input("Enter number of error to correct t: "))
    k = n - r
    ErrorCorrectingCode(r=r, n=n, t=t)
    message = int(input("Enter a message: "), 2)
    error = int(input("Enter an error in binary, 0 if you want to skip: "), 2)
    if error == 0:
        error = get_random_number_of_hamming_weight(n, t)
    config = yaml.safe_load(open('config.yml'))
    codeword, distorted_codeword, error = encode(
        coder_file=config['coder-generator'],
        message=message,
        m_length=k,
        error=error)
    print("Code word: {0:>0{a}b}".format(codeword, a=n))
    print("Distorted code word: {0:>0{a}b}".format(distorted_codeword, a=n))
    print("Error vector: {0:>0{a}b}".format(error, a=n))
    distorted_codeword = int(input("Enter a distorted codeword: "), 2)
    decoded_message = decode(
        parity_check_file=config['decoder-parity-check'],
        n=n,
        syndrome_file=config['decoder-syndrome-decoding'],
        distorted_code=distorted_codeword)
    print("Decoded message: {0:>0{a}b}".format(decoded_message, a=n - r))
