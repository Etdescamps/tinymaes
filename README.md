# TinyMAES

## Description
TinyMAES is a lightweight library for non-linear mono-objective optimisation
written in standard C using only the functions from the standard C library,
so it does not rely on external dependencies.

It uses the MA-ES algorithm that is a variant of [CMA-ES](http://cma.gforge.inria.fr/)
that does not require the eigendecompostion of matrices.

The algorithm MA-ES is described in this article:
> "Simplify Your Covariance Matrix Adaptation Evolution Strategy"
> by Hans-Georg Beyer and Bernhard Sendhoff
> published in IEEE TRANSATIONS ON EVOLUTIONARY COMPUTATION, VOL. 21 no. 5, OCTOBER 2017

## Features
TinyMAES has been written so it can be easily embedded within an application.
It does not use any global variable, so it is possible to run multiple instance
of this algorithm in different threads.

This library contains:
 * A Mersenne Twister Random Number Generator that is based on
   "Tables of 64-bit Mersenne Twisters" by T. Nishimura
   (published on ACM Transactions on Modeling and Computer Simulation 10. (2000))
 * A heap sort function that is especially efficient when lambda >> mu.


