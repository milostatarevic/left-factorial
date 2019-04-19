# Improved Algorithms for Left Factorial Residues

This repository provides the implementation of the algorithms presented in the paper
"Improved Algorithms for Left Factorial Residues" by
Vladica Andrejic, Alin Bostan, and Milos Tatarevic.

## Dependencies

* **NTT**: A library for large integer arithmetic (included in the source)
* [GMP: The GNU Multiple Precision Arithmetic Library](https://gmplib.org)
* [FLINT: Fast Library for Number Theory](https://github.com/wbhart/flint2)
* [NTL: A Library for doing Number Theory](https://www.shoup.net/ntl)

## Usage

### Compute !p mod p for a single prime p

The implementation code is placed in the `single` folder.
After build, you can run, for example:

```
./single/lfr 1099508390819
```

which should return:

```
1099508390819	3851026
```

where `3851026` is `!1099508390819 mod 1099508390819`.

### Compute !p mod p for all primes p in a given range

The implementation code is placed in the `range` folder.
If you want to calculate left factorial residues for all the primes `p <= 1000000007`
by using 6 threads, run:

```
./range/lfrs 0 1000000000 0 6 tree saved backup inverse
```

The results will be saved in `out.remainders.2-1000000007`.
The files stored in the `saved` folder will be used the next time you continue the computation
(starting from 1000000007). If you want to extend the computation up to the first prime greater than
2000000000, run:

```
./range/lfrs 1000000008 2000000000 1000000007 6 tree saved backup inverse
```

The results will be saved in `out.remainders.1000000009-2000000011`.

Note that covering the ranges beyond 2^40 will require disk space measured in tens of terabytes.
In that case, it might be useful to keep `tree`, `saved` and `backup` folders on the separate HDDs.
To reach better performance, it's recommended to keep the `inverse` folder on a SSD.

Note that some values in `defs.h`, such as `MUL_STEP_DEPTH` and `STEP_BITS`,
are calibrated for the setup and the search range we covered (2^34, 2^40).
There are no guarantees these values will provide optimal results for an arbitrary search
range. Always ensure the program is not reporting any errors during the execution.
