# Author: Hubert Kario, Stefan Dordevic
# Released under Gnu GPL v2.1, see LICENSE file for details

try:
    from ecdsa.ellipticcurve import PointJacobi
    JACOBI = True
except ImportError:
    from ecdsa.ellipticcurve import Point
    JACOBI = False
from ecdsa.ellipticcurve import CurveFp, INFINITY
from ecdsa.numbertheory import gcd, jacobi

oneseventwoeight = 1728


def lucasPQ(p, q, n, m):
    """
    nth element of lucas sequence with
    parameters p and q (mod m)
    """
    def half(x):
        if x % 2 == 1: x = x + m
        return (x // 2) % m
    un, vn, qn = 1, p, q
    u = 0 if n % 2 == 0 else 1
    v = 2 if n % 2 == 0 else p
    k = 1 if n % 2 == 0 else q
    n, d = n // 2, p * p - 4 * q
    while n > 0:
        u2 = (un * vn) % m
        v2 = (vn * vn - 2 * qn) % m
        q2 = (qn * qn) % m
        n2 = n // 2
        if n % 2 == 1:
            uu = half(u * v2 + u2 * v)
            vv = half(v * v2 + d * u * u2)
            u, v, k = uu, vv, k * q2
        un, vn, qn, n = u2, v2, q2, n2
    return u, v, k


def lucas_lehmer_riesel_test(n, s, q):
    """
    The "N+1 Test".

    While it builds on top of Lucal Lehmer Riesel test, it's not really
    typical LLR.
    Described in Brillhart, Lehmer, Selfridge (1975) "New Primality Criteria
    and Factorizations of $2^m \pm 1", Theorem 15
    """
    p = q % 2 + 1
    r, rem = divmod(n + 1, s)
    if rem != 0:
        raise ValueError("Invalid certificate")

    if gcd(p, q) != 1:
        raise ValueError("Invalid certificate")
    if jacobi(p ** 2 - 4 * q, n) != -1:
        raise ValueError("Invalid certificate")
    u, v, k = lucasPQ(p, q, (n+1)//2, n)
    if not v == 0:
        raise ValueError("Invalid certificate")
    u, v, k = lucasPQ(p, q, s//2, n)
    if not v != 0:
        raise ValueError("Invalid certificate")
    return r


def pocklington_test(n, s, b):
    """
    The "N-1 Test".

    Described in Brillhart, Lehmer, Selfridge (1975) "New Primality Criteria
    and Factorizations of $2^m \pm 1", Theorem 4
    """
    r, rem = divmod(n - 1, s)
    if rem != 0:
        raise ValueError("Invalid certificate")
    if not s < r:
        raise ValueError("Invalid certificate")
    if pow(b, n-1, n) != 1:
        raise ValueError("Invalid certificate")
    if gcd((pow(b, s, n) - 1) % n, n) != 1:
        raise ValueError("Invalid certificate")
    return r


def curve_test(n, s, w, t, a=None, b=None, j=None):
    """
    Atkin-Goldwasser-Kilian-Morain test.

    the "Elliptic Curve Test"
    """
    if a is None:
        a = (3 * j * (oneseventwoeight - j)) % n
    if b is None:
        b = (2 * j * (oneseventwoeight - j)**2) % n

    r, rem = divmod(n + 1 - w, s)
    if rem != 0:
        raise ValueError("Invalid certificate")
    l = (t ** 3 + a * t + b) % n
    a = (a * l ** 2) % n
    b = (b * l ** 3) % n
    x = (t * l) % n
    y = (l ** 2) % n

    curve = CurveFp(n, a, b)
    if JACOBI:
        p1 = PointJacobi(curve, x, y, 1)
    else:
        p1 = Point(curve, x, y)
    p2 = p1 * s
    if p2 == INFINITY:
        raise ValueError("Invalid certificate")
    p3 = p2 * r
    if p3 != INFINITY:
        raise ValueError("Invalid certificate")
    return r
