# Refined Humbert Invariants

This repository contains Julia code for computing refined Humbert invariants for a random set of polarizations given a prime $p\equiv 11\mod 12$. 

## Code Organization

### Genus of the canonical RHI

For prime $p\equiv 11 \mod 12$, the canonical refined Humbert invariant is the one corresponding to the identity polarization 

$$q_{(E \times E, \theta)}(D)=t_0^2+4\left(t_1^2 + t_1t_3 + t_2^2 + t_2t_4 + \frac{p+1}{4}t_3^2 + \frac{p+1}{4}t_4^2\right).$$

The [`genus`](./genus/) directory contains the code and data about the genus of this 5-ary form for primes of interest.

### Random polarizations

We adapt the `RandomPolarisation` function of [KLPT2](https://github.com/KLPT2/KLPT2) to compute a set of random polarizations on a fixed superspecial surface $$E \times E$$.

The [`polz`](./polz/) directory contains the code and data about such polarizations for primes of interest. 

### Degree maps

The refined Humbert invariants which represent 1 can be expressed as

$$ q_{E\times E, \theta}(D) = X^2 + 4\cdot \deg$$

where $\deg$ is a 4-ary quadratic form.

The [`deg`](./deg/) directory contains the code and data about degree maps and $\max(\min(\deg))$ values for the primes of interest for random polarizations.

### Automorphism groups

The refined Humbert invariants which do not represent 1, can be used to infer distribution of automorphism group of genus 2 curve $C$ because

$$i(\mathrm{Aut}(C)) - 1 = r_4(q_C)$$

where $i(\mathrm{Aut}(C))$ is the number of involutions of the automorphism group of $C$, and $r_4(q_C)$ is the number of ways the irreducible refined Humbert invariant $q_C$ represents 4.

The [`aut`](./aut/) directory contains the code and data about the automorphism group distribution for the primes of interest.

### Splittings

The [`split`](./split/) directory contains the data about the splittings.

### Superspecial polarizations

This algorithm classifies principal polarizations by enumerating unimodular positive definite binary Hermitian forms over the fixed maximal order $O \subset \mathrm{B}_p$. We adapt lattice techniques to switch from binary Hermitian module to integral lattice of rank 8.

Integral lattices come equipped with 4 trace forms in 8 variables (auxiliary forms) that are respecting to both the quadratic form and the O-structure. We reduce the auxiliary forms by LLL, and then applying exact simultaneous Plesken–Souvignier isometry. The enumeration bound is increased dynamically, and the computation terminates once the number of non-isomorphic principal polarizations reaches the class number $H(p)$. 

We adapt the ideas given in the research paper by [Greenberg and Voight](http://arxiv.org/abs/1209.2460v1), and [Sarah Chisholm's dissertation](http://hdl.handle.net/11023/1920).

The [`optimal`](./optimal/) directory contains code and data about polarizations and corresponding refined Humbert invariants.

## System Requirements

The code was written and tested on a system with the following specifications:

CPU: 12th Gen Intel i7-12800HX (24)

Memory: 15840MiB

OS: Ubuntu 22.04.5 LTS on Windows 11 x86_64 (WSL2)

Software: Julia Version 1.12.4 with Oscar Version 1.6.0
