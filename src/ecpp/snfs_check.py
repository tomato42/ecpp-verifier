import math


def snfs_vulnerable(n, distance=2**64):
    """
    Check if number `n` is further away from any exact power than `distance`.

    Test to verify if the number is vulnerable to Special Number Field Sieve
    attack on discrete logarithm problem.

    Return True if it is vulnerable, False if it is not vulnerable.
    """
    max_base = 100
    # For bases smaller than max_base, just calculate the value by calculating
    # logarithm
    for base in range(2, max_base):
        log = math.log(n, base)
        if min(abs(base**int(math.floor(log)) - n),
               abs(base**int(math.ceil(log)) - n)) <= distance:
            return True

    # for higher bases, the result of logarithm
    # (the exponent we're searching for)
    # changes only a little for every base increase,
    # and we are only interested in integer bases and exponents,
    # so reverse the search -
    # look for bases for which the integer exponent will lie near n
    max_exponent = int(math.log(n, max_base - 1)) + 1
    for exponent in range(max_exponent, 2, -1):
        # e = log(n, b) = log(n, 2) / log(b, 2)
        # log(b, 2) = log(n, 2) / e
        # b = 2**(log(n, 2) / e)
        base_high = 2**int(math.ceil(float(n.bit_length()) / exponent))
        base_low = 2**int(math.floor((float(n.bit_length()) - 1) / exponent))
        # assert base_low**exponent <= n <= base_high**exponent
        while base_high != base_low + 1:
            attempt_at_base = (base_low + base_high) // 2
            attempt_at_n = attempt_at_base**exponent
            if attempt_at_n > n:
                base_high = attempt_at_base
            elif attempt_at_n < n:
                base_low = attempt_at_base
            else:
               return True  # base**e = n = attempt_at_n
        if min(abs(base_high**exponent - n),
               abs(base_low**exponent - n)) <= distance:
            return True
    return False


assert snfs_vulnerable((2**8192) + 1)
assert snfs_vulnerable((2**8192) - 1)
assert not snfs_vulnerable((2**8192) - 2**64 - 1)
assert snfs_vulnerable((2**8192) + 2**64)
assert not snfs_vulnerable(801**921 + 2**64 + 1)
assert snfs_vulnerable(801**921 + 2**64)
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
