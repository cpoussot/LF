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

The three following exercises are set up to discover the LF and some of its features. It can be used to discover step by step how to set up and solve the interpolation problem.

## Exercise \#1: first try with Loewner

Define two toy rational transfer functions $G_1(s) = \frac{1}{2+s}$ and  $G_2(s) = \frac{s+1}{s+5}$ as a Matlab handle functions.
```Matlab
G1 = @(s) 1/(s+2);
G2 = @(s) (s+1)/(s+5);
G  = G1; % or G = G2;
nu = 1;
ny = 1;
```

### Construct data (interpolation points and its evaluation)

Now, create 20 interpolation points, and split them in $\lambda$ and $\mu$, e.g.
```Matlab
si  = linspace(1e-2,100,20);
la_ = si(1:2:end);
mu_ = si(2:2:end);
```

Now construct the data, being the evaluation of $G$ at interpolation points $\lambda$'s and $\mu$'s (for tangential direction being just ones, as the system is SISO). 
```Matlab
k   = length(la_);
q   = length(mu_);
R   = ones(nu,k);
L   = ones(q,ny);
for ii = 1:k; W(1:ny,1:nu,ii) = G(la_(ii)); end
for ii = 1:q; V(1:ny,1:nu,ii) = G(mu_(ii)); end
```

### Construct the rational interpolant (apply Loewner tangential)

Now you may run the tangential Loewner approximation with an SVD threshold quite low, e.g. a rolerance of $10^{-12}$.
```Matlab
opt         = [];
opt.target  = 1e-12;
[htng,itng] = lf.loewner_tng(la_,mu_,W,V,R,L,opt);
```

Examinate the outputs `htng` and `itng`, and conclude for either $G_1$ and $G_2$,
- on the approximation order $r$ and McMillan degree $\nu$;
- on the accuracy of the eigenvalues of the obtained model;
- on the normalized singular values of the Loewner pencil;
- on frequency gain and phase response (by plotting Bode and/or Nyquist);
- on the Sylvseter equations property;
- on the informations contained in `itng` (have a deeper look at these informations).

## Exercise \#2: dig a bit more...

### Construct data

Now, create interpolation points along the imaginary axis, and split them in $\lambda$ and $\mu$ (and close them by conjugation), e.g.
```Matlab
nip = 10;
la_ = (logspace(-2,2,nip))*1i;    
la_ = sort([la_ conj(la_)]);
mu_ = (logspace(-2,2,nip)+.1)*1i; 
mu_ = sort([mu_ conj(mu_)]);
```

Construct the new data `W` and `V` just as in Exercise \#1.

### Construct the rational interpolant (apply Loewner tangential)

Now compute the interpolant without enforcing model realness.
```Matlab
opt         			= [];
opt.target  			= 1e-12;
opt.real    			= false;
[htng_cplx,itng_cplx] 	= lf.loewner_tng(la_,mu_,W,V,R,L,opt);
```

Examinate the outputs `htng_cplx` and `itng_cplx`, and conclude for either $G_1$ and $G_2$,
- on the approximation order $r$ and McMillan degree $\nu$;
- on the accuracy of the eigenvalues of the obtained model;
- on the normalized singular values of the Loewner pencil;
- on frequency gain and phase response (by plotting Bode and/or Nyquist);
- on the Sylvseter equations property;
- on the informations contained in `itng_cplx` (have a deeper look at these informations).

Notice that the matrices and values are all complex.

### Construct the real rational interpolant (apply Loewner tangential)

Now compute the interpolant by enforcing model realness.
```Matlab
opt         			= [];
opt.target  			= 1e-12;
opt.real    			= true;
[htng_real,itng_real] 	= lf.loewner_tng(la_,mu_,W,V,R,L,opt);
```

Examinate the outputs `htng_real` and `itng_real`, and conclude
- on the approximation order $r$ and McMillan degree $\nu$;
- on the accuracy of the eigenvalues of the obtained model;
- on the normalized singular values of the Loewner pencil;
- on the Sylvseter equations property;
- on frequency gain and phase response (by plotting Bode and/or Nyquist);
- on the informations contained in `itng_real` (have a deeper look at these informations).

Notice that the matrices and values are now real.


## Exercise \#3: about reduction, passive models and pH structures

### Run the demo file `demo_ph.m`

Now we suggest to try the script `demo_ph.m`. The latter uses the very same LF but with an additional step to ensure passivity and port Hamiltonian structre. Please spend some time to test the different models proposed
```Matlab
CAS = 'siso_passive_simple';   r = 1e-12
CAS = 'siso_passive_aca';      r = 1e-12
CAS = 'siso_passive_gugercin'; r = 5
CAS = 'mimo_passive_aca';      r = 5
```

Examinate the outputs of the passive rational approximant `hloep`, `info_loep` and `info_loeph`, and conclude
- on the approximation order $r$ and McMillan degree $\nu$;
- on the accuracy of the eigenvalues of the obtained model;
- on the normalized singular values of the Loewner pencil;
- on the Sylvseter equations property;
- on frequency gain and phase response (by plotting Bode and/or Nyquist);
- on the spectral zeros, positive spectral zeros;
- on the Cholesky projector;
- on the informations contained in `info_loep` (have a deeper look at these informations).
- on the informations contained in `info_loeph` (have a deeper look at these informations, especially on the pH matrices)


## Exercise \#4: simple mass-spring-damper example

Now use the example suggested in the Part 3 of the slides. The MSD is given the following equation
$$
\begin{array}{rcl}
E\dot x(t) &=& Ax(t)+Bu(t)\\
y(t) &=& Cx(t)
\end{array}
$$


### Define a state-space model of the MSD

```Matlab
% MSD parameters
m 	= 10; % mass
k 	= 1;  % spring
b 	= 3;  % damper
% System properties
n   = 2;
ny  = 1;
nu  = ny;
% State-space X=[x dx]
E   = eye(n);
A   = [0 1; -k/m -b/m];
B   = [0; 1/m];
C   = [0 1];
D   = 0;
Hss = dss(A,B,C,D,E);
hss = @(s) C*((s*E-A)\B) + D;
```

### Apply LF 


### Go to the time-domain simulations




## Exercise \#5: go far away and measure the flexibility of LF 

- Play with target order.
- Change the original function and try replace with an irrational one such as $G(s)=\frac{1}{se^{-s}+1}$.



