# Author: Hubert Kario, Stefan Dordevic
# Released under Gnu GPL v2.1, see LICENSE file for details

from ecdsa.ellipticcurve import CurveFp, Point, INFINITY
from ecdsa.numbertheory import gcd, jacobi

oneseventwoeight = 1728

def fibonacci(n, mod):
#    print("calculating for: " + str(n))
    if n == 0:
        return 0
    # en.wikipedia.org/wiki/Fibonacci_number#Matrix_form
    v1, v2, v3 = 1, 1, 0 # init matrix [[1, 1], [1, 0]]
    for i, rec in enumerate(bin(n)[3:]):  # perform fast exponentiation of the matrix (quickly raise it to the nth power)
#        print(i)
        calc = v2*v2 % mod
        v1, v2, v3 = (v1*v1+calc) % mod, ((v1+v3)*v2) % mod, (calc+v3*v3) % mod
        if rec=='1':
            v1, v2, v3 = (v1+v2) % mod, v1, v2
    return v2

fib_matrix = [[1,1],
              [1,0]]

def matrix_square(A, mod):
    return mat_mult(A,A,mod)


def mat_mult(A,B, mod):
  if mod is not None:
    return [[(A[0][0]*B[0][0] + A[0][1]*B[1][0])%mod, (A[0][0]*B[0][1] + A[0][1]*B[1][1])%mod],
            [(A[1][0]*B[0][0] + A[1][1]*B[1][0])%mod, (A[1][0]*B[0][1] + A[1][1]*B[1][1])%mod]]


def matrix_pow(M, power, mod):
    #Special definition for power=0:
    if power <= 0:
      return M

    powers =  list(reversed([True if i=="1" else False for i in bin(power)[2:]])) #Order is 1,2,4,8,16,...

    matrices = [None for _ in powers]
    matrices[0] = M

    for i in range(1,len(powers)):
        matrices[i] = matrix_square(matrices[i-1], mod)


    result = None

    for matrix, power in zip(matrices, powers):
        if power:
            if result is None:
                result = matrix
            else:
                result = mat_mult(result, matrix, mod)

    return result


def fibonacci2(n, mod):
    return matrix_pow(fib_matrix, n, mod)[0][1]


assert fibonacci(0, 2**1024) == 0
assert fibonacci(1, 2**1024) == 1
assert fibonacci(2, 2**1024) == 1
assert fibonacci(3, 2**1024) == 2
assert fibonacci(4, 2**1024) == 3
assert fibonacci(5, 2**1024) == 5
assert fibonacci(12, 2**1024) == 144
assert fibonacci(400, 2**1024) == 176023680645013966468226945392411250770384383304492191886725992896575345044216019675
assert fibonacci(1000, 2**1024) == 43466557686937456435688527675040625802564660517371780402481729089536555417949051890403879840079255169295922593080322634775209689623239873322471161642996440906533187938298969649928516003704476137795166849228875
assert fibonacci(12, 100) == 144 % 100
assert fibonacci(400, 100) == 176023680645013966468226945392411250770384383304492191886725992896575345044216019675 % 100
assert fibonacci(1000, 100) == 43466557686937456435688527675040625802564660517371780402481729089536555417949051890403879840079255169295922593080322634775209689623239873322471161642996440906533187938298969649928516003704476137795166849228875 % 100
assert fibonacci(12, 11) == 144 % 11
assert fibonacci(400, 11) == 176023680645013966468226945392411250770384383304492191886725992896575345044216019675 % 11
assert fibonacci(1000, 11) == 43466557686937456435688527675040625802564660517371780402481729089536555417949051890403879840079255169295922593080322634775209689623239873322471161642996440906533187938298969649928516003704476137795166849228875 % 11
assert fibonacci(10**10, 87) == fibonacci2(10**10, 87)
assert fibonacci(10**11, 83) == fibonacci2(10**11, 83)
assert fibonacci(10**19, 1000000007) == fibonacci2(10**19, 1000000007)

def lucas(i, mod):
    return (fibonacci(i - 1, mod) + fibonacci(i + 1, mod)) % mod

def lucas2(i, mod):
    return (fibonacci2(i - 1, mod) + fibonacci2(i + 1, mod)) % mod

#assert lucas(0) == 2
assert lucas(1, 2**1024) == 1
assert lucas(2, 2**1024) == 3
assert lucas(10, 2**1024) == 123
assert lucas(10, 12) == 123 % 12
assert lucas(1000, 12) == lucas(1000, 2**1024) % 12
assert lucas(10, 83) == 123 % 83
assert lucas(1000, 83) == lucas(1000, 2**1024) % 83
assert lucas(10**19, 1000000005) == lucas2(10**19, 1000000005)

def lucasPQ(p, q, n, m):
#    print(p, q, n)
    # nth element of lucas sequence with
    # parameters p and q (mod m)
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
    # the "N+1 Test"
    # while it builds on top of Lucal Lehmer Riesel test, it's not really
    # typical LLR
    # Described in Brillhart, Lehmer, Selfridge (1975) "New Primality Criteria
    # and Factorizations of $2^m \pm 1", Theorem 15
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
    # the "N-1 Test"
    # Described in Brillhart, Lehmer, Selfridge (1975) "New Primality Criteria
    # and Factorizations of $2^m \pm 1", Theorem 4
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
    # Atkin-Goldwasser-Kilian-Morain test
    # the "Elliptic Curve Test"
    
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
    p1 = Point(curve, x, y)
    p2 = p1 * s
    if p2 == INFINITY:
        raise ValueError("Invalid certificate")
    p3 = p2 * r
    if p3 != INFINITY:
        raise ValueError("Invalid certificate")
    return r
