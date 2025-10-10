# Loewner Framework (LF)

## Overview

The Loewner Framework (LF) is a data-driven model identification and reduction technique that takes its root in ["Rational interpolation and state-variable realizations"](https://doi.org/10.1016/0024-3795(90)90140-8), by B.D.O. Anderson and A.C. Antoulas, where connection between the Loewner matrix and rational barycentric form was highlighted. Later in ["A framework for the solution of the generalized realization problem"](https://doi.org/10.1016/j.laa.2007.03.008), by A.J. Mayo and A.C. Antoulas, authors extented this framework to connect the Loewner (and shifted Loewner) matrix to the property of minimal realization, bridging the gap between rational interpolation, realization and system theory. 

The Loewner framework (LF) can thus be viewed as a data-driven model identification, reduction and approximation tool. It provides a solution through rational approximation by means of interpolation. Using only measured data, the LF constructs surrogate models directly and with low computational effort (e.g. using only basic linear algebra operations and SVD-like deocompositions). For recent tutorial papers on LF applied to linear dynamical systems, we refer the reader to ["A Tutorial Introduction to the Loewner Framework for Model Reduction"](https://doi.org/10.1137/1.9781611974829.ch8), by A.C. Antoulas, A.C. Ionita and S. Lefteriu.

It is suited to approximate, from any data-set:
- an interpolatory model given with a realization,
- with minimal complexity.

## Contribution (dis)claim

The present page does not aim neither at describing this framework in details nor providing a novel contribution, but rather to provide elementary ingredients to test and teach the LF thanks to MATLAB codes and the package `+lf`. The latter includes the following features and functions:
- approximate any SISO to MIMO data-set with a rational functions;
- with controlled complexity; 
- with passivity preservation;
- and with port Hamiltonian realization.

The material contained in this package is mostly used to feel the LF and for education purpose. Additional features to deal with passivity and port-Hamiltonian  systems are also given.

## Main reference (from the author)

```
@article{GPVA:2022,
	Author	= {I.V. Gosea and C. Poussot-Vassal and A.C. Antoulas},
	Doi 	= {10.1016/bs.hna.2021.12.015},
	Journal = {Handbook in Numerical Analysis, Vol. 23, January 2022, pp. 499-530},
	Number 	= {},
	Pages 	= {499-530},
	Title 	= {{Data-driven modeling and control of large-scale dynamical systems in the Loewner framework}},
	Volume 	= {23},
  	Month   = {January},
	Year 	= {2022},
	Note    = {\url{https://doi.org/10.1016/bs.hna.2021.12.015}}, 
}
```

# The "LF" MATLAB package 

The code (`+lf` folder)  provided in this GitHub page is given for open science purpose. Its principal objective is to accompany the the reasers, and thus aims at being as educative as possible rather than industry-oriented. Evolutions (numerical improvements) may come with time. Please, cite the reference above if used in your work and do not hesitate to contact us in case of bug of problem when using it. Below we present an example of use, then functions list are given.

## A simple MATLAB code example

Here is a simple code that describes how to deploy the LF. Code below is `demo.m`.

First add the path where the `+mlf` package is located.

```Matlab
%addpath("location_of_lf") % Add the location of the +lf package
```

## Functions description

Soon


## Feedbacks

Please send any comment to C. Poussot-Vassal (charles.poussot-vassal@onera.fr) if you want to report any bug or user experience issues.

## Disclaimer

Once again, this deposit consitutes a research code that accompany the paper mentionned above. It is not aimed to be included in any third party software without the consent of the authors. Authors decline responsabilities in case of problem when applying the code.

Notice also that pathological cases may appear. A more advanced code, to deal with practical and theoretical issues/limitations is currently under developpement by the authors.


# Exercises suggestions

## Exercise \#1 

Define the rational transfer function as a handle function 
$$
G(s) = \frac{s+1}{5+s}
$$
as 

```Matlab
G = @(s) (s+1)/(s+5)
```

