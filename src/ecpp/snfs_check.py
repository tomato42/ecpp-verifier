import math
import gmpy2


def snfs_vulnerable(n, distance=2**64):
    """
    Check if number `n` is further away from any exact power than `distance`.

    Test to verify if the number is vulnerable to Special Number Field Sieve
    attack on discrete logarithm problem.

    Return True if it is vulnerable, False if it is not vulnerable.
    """
    gmpy2.get_context().precision = 1024

    max_base = 100
    # For bases smaller than max_base, just calculate the value by calculating
    # logarithm
    for base in range(2, max_base):
        exponent = int(math.log(n, base))
        if abs(pow(base, exponent) - n) <= distance:
            return True
        # int() rounds down, try value rounded up
        if abs(pow(base, exponent + 1) - n) <= distance:
            return True
    # for higher bases, the result of logarithm (the exponent we're
    # searching for) changes only a little for every base increase, and we
    # are only interested
    # in integer bases and exponents, so reverse the search - look for bases
    # for which the integer exponent will lie near n
    max_exponent = int(math.log(n, max_base - 1)) + 1
    for exponent in range(max_exponent, 2, -1):
        # compute n ^ (1/exponent), but since n is huge, it causes float
        # overflow so perform the same calculation but through natural logs
        base = int(gmpy2.exp(gmpy2.log(n) / gmpy2.mpfr(exponent)))
        if abs(pow(base, exponent) - n) <= distance:
            return True
        # int() rounds down, so try the value rounded up
        if abs(pow(base + 1, exponent) - n) <= distance:
            return True

    return False


assert snfs_vulnerable((2**8192) + 1)
assert snfs_vulnerable((2**8192) - 1)
assert not snfs_vulnerable((2**8192) - 2**64 - 1)
assert snfs_vulnerable(2**1024 - 2**63)
assert snfs_vulnerable(3**3072)
assert snfs_vulnerable(101**100)
assert not snfs_vulnerable(800**300 + 2**64 + 1)
assert snfs_vulnerable(800**300 + 2**64)
assert not snfs_vulnerable(800**300 - 2**64 - 1)
assert snfs_vulnerable(800**300 - 2**64)
assert not snfs_vulnerable(2**1024 + 2**64 + 1)
assert snfs_vulnerable(2**1024 + 2**64)
assert not snfs_vulnerable(2**1024 - 2**64 - 1)
assert snfs_vulnerable(2**1024 - 2**64)
# also check with an actual safe-prime:
assert snfs_vulnerable(2**1024 - 1093337)
