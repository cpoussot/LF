# /!\ UNDER CONSTRUCTION

# Loewner Framework (LF)

### Overview

The Loewner Framework (LF) is a data-driven model identification and reduction technique that takes its root in ["Rational interpolation and state-variable realizations"](https://doi.org/10.1016/0024-3795(90)90140-8), by B.D.O. Anderson and A.C. Antoulas, where connection between the Loewner matrix and rational barycentric form was highlighted. Later in ["A framework for the solution of the generalized realization problem"](https://doi.org/10.1016/j.laa.2007.03.008), by A.J. Mayo and A.C. Antoulas, authors extented this framework to connect the Loewner (and shifted Loewner) matrix to the property of minimal realization, bridging the gap between rational interpolation, realization and system theory. 

The Loewner framework (LF) can thus be viewed as a data-driven model identification, reduction and approximation tool. It provides a solution through rational approximation by means of interpolation. Using only measured data, the LF constructs surrogate models directly and with low computational effort (e.g. using only basic linear algebra operations and SVD-like deocompositions). For recent tutorial papers on LF applied to linear dynamical systems, we refer the reader to ["A Tutorial Introduction to the Loewner Framework for Model Reduction"](https://doi.org/10.1137/1.9781611974829.ch8), by A.C. Antoulas, A.C. Ionita and S. Lefteriu.

It is suited to approximate, from any data-set:
- an interpolatory model given with a realization,
- with minimal complexity.

### Package purpose

The present page does not aim at describing this framework in details, but rather to provide elementary ingredients to test and teach the LF thanks to Matlab code. The latter includes the following features:
- to approximate any SISO to MIMO data-set with a rational functions;
- with controlled complexity; 
- and with passivity preservation;

### Functions descriptions





