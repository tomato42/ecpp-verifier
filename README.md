Collection of primality certificates and an open source tool
that can be used to verify them.

Application can verify [Atkin-Goldwasser-Kilian-Morain Certificate](http://mathworld.wolfram.com/Atkin-Goldwasser-Kilian-MorainCertificate.html),
[Pocklington certificate](http://mathworld.wolfram.com/PocklingtonsTheorem.html)
and
[Brillhart, Lehmer, Selfridge certificate](http://www.ams.org/mcom/1975-29-130/S0025-5718-1975-0384673-1/S0025-5718-1975-0384673-1.pdf) (Theorem 15) based
primality proofs, commonly generated by the [Primo](http://www.ellipsa.eu/public/primo/primo.html)
application.

Read more on [primality testing](http://cr.yp.to/primetests.html)

When verifying primality certificates, it also checks if the number is not
vulnerable to
[Special Number Field Sieve](https://en.wikipedia.org/wiki/Special_number_field_sieve),
as primes of this form have effective bit size reduced by a third
(i.e. 3072 bit SNFS-vulnerable prime is about as easy to attack as a 2048 bit
prime).

# Installation

## From distribution files

Download the most recent release `.whl` file from github:

    curl -s https://api.github.com/repos/tomato42/ecpp-verifier/releases/latest \
    | grep "browser_download_url.*whl" \
    | cut -d : -f 2,3 \
    | tr -d \" \
    | wget -i -

Or manually, by visiting
[latest release](https://github.com/tomato42/ecpp-verifier/releases/latest)
and downloading the `.whl` file from there.

Install it using `pip`:

    pip3 install ecpp-*.whl

Verify that it can be executed:

    ecpp --help

Install `gmpy2` package to double the performance of certificate verification:

    pip3 install gmpy2

## From sources

Clone this repository:

    git clone https://github.com/tomato42/ecpp-verifier.git
    cd ecpp-verifier

Install dependencies (for example on Fedora):

    dnf install python3-ecdsa python3-gmpy2

Or from PyPI:

    pip3 install ecdsa[gmpy2]

(Note: as `gmpy2` is a binary package you will need to install development
headers for python and the gmp library. Alternatively, you can skip
installation of `gmpy2` at the cost of halved performance).

Make sure you have installed `ecdsa` package version 0.15 or newer. Older
versions have significant performance issues and certificate verification
will take ages.

Run `ecpp` for the first time:

    PYTHONPATH=src ./ecpp --help

# Usage

## Matching certificates to primes in OpenSSH moduli file

To check if you have primality certificates for all the primes in your OpenSSH
moduli file, you can use ecpp with `-m` switch:

    ecpp -m /etc/ssh/moduli

This will succeed for example for moduli file released with OpenSSH 8.2p1,
listing certificates for each prime.

## Verifying primality certificates for primes in moduli file

To verify the matching certificates you can combine the `-m` switch with
the `-v` switch:

    ecpp -m /etc/ssh/moduli -v

This will succeed if the script finds matching certificates and verifies
them as valid.

Note: it will require significantly more time to execute than just the `-m`
option. It's also a single-threaded process, see
[#12](https://github.com/tomato42/ecpp-verifier/issues/12).

## Generating primality certificates for primes in moduli file for OpenSSH

If there are some primes without primality certificates, you can generate
input files for Primo into `in/` directory.

    ecpp -m /etc/ssh/moduli -p

Now, open Primo downloaded from link above, extract archive and start GUI on
a reasonably powerful machine (at this moment, Primo can work with up to
64 cores).

 * From *Menu*, select *Setup...*, set number of cores your system have
   (hyper-threading is not much useful)
 * In tab *Certification*, select *3000 dd* (decimal digits) in *Trial
   Division Parameters* and click *Build prime table* button.
 * Then click *Load*. Select all the `.in` files in
   the `in/` directory created by the previous step and click *Open*.
 * After long time, you will get certificates in `*.out` files in the same
   directory.

## Verify primality certificates

In previous step, we got certificates for primes. Now we need to verify them.
This can be done with the following command for one certificate:

    ecpp -i in/primo-B412D0397A9D9-07E.out

The job can be simply parallelized so if we want to verify all
the certificates we got, we can use GNU parallel to get results in parallel,
in this example using 16 parallel processes:

    parallel -j16 --halt now,fail=1 "echo {} && ecpp -i {}" ::: in/*.out

or with progress reporting:

    parallel --eta --progress --bar -j16 --halt now,fail=1 "echo {} && ecpp -i {}" ::: in/*.out

(Change the parameter to the `-j` option if you have fewer or more than 16
cores)

Now, we can add the primality certificates to `src/ecpp/certificates/`
directory.

## Matching certificates to primes in OpenSSH moduli file again

Running `ecpp` again as in the first example should confirm we have
a certificate for each prime in the moduli file now.

    ecpp -m /etc/ssh/moduli
