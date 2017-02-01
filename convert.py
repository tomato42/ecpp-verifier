# Author: Hubert Kario, Stefan Dordevic
# Released under Gnu GPL v2.1, see LICENSE file for details

from __future__ import division
import sys
import hashlib

"""
Functions for hashing the prime number  

"""

def numBits(n):
    """Return number of bits necessary to represent the integer in binary"""
    if n == 0:
        return 0
    if sys.version_info < (2, 7):
        # bit_length() was introduced in 2.7, and it is an order of magnitude
        # faster than the below code
        return len(bin(n))-2
    else:
        return n.bit_length()

def numBytes(n):
    """Return number of bytes necessary to represent the integer in bytes"""
    if n == 0:
        return 0
    bits = numBits(n)
    return (bits + 7) // 8

def numberToByteArray(n, howManyBytes=None, endian="big"):
    """
    Convert an integer into a bytearray, zero-pad to howManyBytes.

    The returned bytearray may be smaller than howManyBytes, but will
    not be larger.  The returned bytearray will contain a big- or little-endian
    encoding of the input integer (n). Big endian encoding is used by default.
    """
    if howManyBytes == None:
        howManyBytes = numBytes(n)
    if endian == "big":
        return bytearray((n >> i) & 0xff
                         for i in reversed(range(0, howManyBytes*8, 8)))
    elif endian == "little":
        return bytearray((n >> i) & 0xff
                         for i in range(0, howManyBytes*8, 8)) 
    else:
        raise ValueError("Only 'big' and 'little' endian supported")

def sha256_of_number(int_prime_moduli):
    """use sha256 to hash the prime

    [Candidate]
    N=<hex>
    """
    n = int_prime_moduli
    n = numberToByteArray(n)
    n = hashlib.sha256(n).hexdigest()
    return n
