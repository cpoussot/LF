# Loewner Framework (LF)

## Overview

The Loewner Framework (LF) is a data-driven model identification and reduction technique that takes its root in ["Rational interpolation and state-variable realizations"](https://doi.org/10.1016/0024-3795(90)90140-8), by B.D.O. Anderson and A.C. Antoulas, where connection between the Loewner matrix and rational barycentric form was highlighted. Later in ["A framework for the solution of the generalized realization problem"](https://doi.org/10.1016/j.laa.2007.03.008), by A.J. Mayo and A.C. Antoulas, authors extented this framework to connect the Loewner (and shifted Loewner) matrix to the property of minimal realization, bridging the gap between rational interpolation, realization and system theory. 

The Loewner framework (LF) can thus be viewed as a data-driven model identification, reduction and approximation tool. It provides a solution through rational approximation by means of interpolation. Using only measured data, the LF constructs surrogate models directly (i.e. without iterations) and with low computational effort (e.g. using only basic linear algebra operations and SVD-like deocompositions). For recent tutorial papers on LF applied to linear dynamical systems, we refer the reader to ["A Tutorial Introduction to the Loewner Framework for Model Reduction"](https://doi.org/10.1137/1.9781611974829.ch8), by A.C. Antoulas, A.C. Ionita and S. Lefteriu. Last, we refer to ["Identification of port-Hamiltonian systems from frequency response data"](https://doi.org/10.1016/j.sysconle.2020.104741), by P. Benner, P. Goyal and P. Van-Dooren for an extention to passive and pH-structures.


It is suited to approximate, from any data-set:
- an interpolatory model given with a realization,
- with minimal complexity,
- and with potentially passivity and pH forlm preservation.

## Contribution (dis)claim

The present page does not aim neither at describing this framework in details nor providing a novel contribution, but rather to provide elementary ingredients to test and teach the LF thanks to MATLAB codes and the suggested package `+lf`. The latter includes the following principal features and functions:
- approximate any SISO to MIMO data-set with a rational functions;
- with controlled complexity; 
- with passivity preservation;
- and with port Hamiltonian realization.

The material contained in this package is mostly used to feel the LF and for education purpose. Additional features to deal with passivity and port-Hamiltonian  systems are also given.

## Main reference (from the author)

In addition to the highly relevant references mentionned in the introduction, interested reader may also refer to the paper ["Data-driven modeling and control of large-scale dynamical systems in the Loewner framework"](https://doi.org/10.1016/bs.hna.2021.12.015), by I.V. Gosea, C. Poussot-Vassal and A.C. Antoulas, which details some recent advances.

```
@article{GPVA:2022,
	Author	= {I.V. Gosea and C. Poussot-Vassal and A.C. Antoulas},
	Doi 	= {10.1016/bs.hna.2021.12.015},
	Journal = {Handbook in Numerical Analysis},
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

The code (`+lf` folder)  provided in this GitHub page is given for open science purpose. Its principal objective is to accompany the readers, and thus aims at being as educative as possible rather than industry-oriented. Evolutions (numerical improvements) may come with time. Please, cite the reference above if used in your work and do not hesitate to contact us in case of bug of problem when using it. Below we present an example of use, then functions list are given.

Moreover, for more numerically robust and involved implementation and features, we invite reader and users to refer to the [MDSPACK](https://mordigitalsystems.fr/static/mdspack_html/MDSpack-guide.html) library by [MOR Digital Systems](https://mordigitalsystems.fr).

## A simple MATLAB code example

Here is a simple code that describes how to deploy the LF. Code below is `demo.m`.

First add the path where the `+mlf` package is located.

```Matlab
%addpath("location_of_lf") % Add the location of the +lf package
```

## Functions description

Please check help in functions
```Matlab
help lf.loewnerMatrix
help lf.loewner_tng
help lf.loewner_blk
help lf.loewner_passive
```

## Feedbacks

Please send any comment to C. Poussot-Vassal (charles.poussot-vassal@onera.fr) if you want to report any bug or user experience issues.

## Disclaimer

Once again, this deposit consitutes a research code that accompany the paper mentionned above. It is not aimed to be included in any third party software without the consent of the authors. Authors decline responsabilities in case of problem when applying the code.

Notice also that pathological cases may appear. A more advanced code, to deal with practical and theoretical issues/limitations is currently under developpement by the authors.


# Exercises suggestions

The following exercises are set up to discover the LF and some of its features. It can be used to discover step by step how to set up and solve the interpolation problem.

## Exercise \#1: simple mass-spring-damper example (MSD)

### Introduction

Let us use the simple sinlge mass-spring-damper (MSD) system. The MSD is given the by following state-space equations,
```math
\begin{array}{rcl}
E\dot{\mathbf{x}}(t) &=& A{\mathbf{x}}(t)+Bu(t)\\
y(t) &=& C{\mathbf{x}}(t)
\end{array}
```
where ($x$ and $\dot x$ respectively are the position and the velocity of the mass),
```math
{\mathbf{x}} = \left(\begin{array}{c} x\\ \dot x \end{array}\right)  ,\,
E = I_2 ,\, 
A = \left(\begin{array}{cc} 0 & 1\\ -\frac{k}{m} & -\frac{b}{m} \end{array}\right) ,\,
B = \left(\begin{array}{c} 0\\ \frac{1}{m} \end{array}\right)  ,\,
C = \left(\begin{array}{cc} 0 & 1 \end{array}\right).
```
Notice that this model is passive (but not stricty). Now remark that it admits a port Hamiltonian (pH) form, e.g. with these equations,
```math
\begin{array}{rcl}
M\dot{\mathbf{x}}(t)&=&(J-R)Q{\mathbf{x}}(t)+(G-P)u(t) \\
y(t)&=&(G+P)^\top Q{\mathbf{x}}(t)+(N+S)u(t)
\end{array}
```
where ($\alpha_x$ and $\alpha_v$ are the energy variables), 
```math
{\mathbf{x}} = \left(\begin{array}{c} \alpha_x\\ \alpha_v \end{array}\right)  ,\,
M = I_2 ,\, 
J = \left(\begin{array}{cc} 0 & \frac{1}{m}\\ -\frac{1}{m} & 0 \end{array}\right) ,\,
R = \left(\begin{array}{cc} 0 & 0\\ 0 & \frac{b}{m^2} \end{array}\right) ,\,
Q = \left(\begin{array}{cc} k & 0\\ 0 & m \end{array}\right) ,\,
G = \left(\begin{array}{c} 0\\ \frac{1}{m} \end{array}\right) ,\,
P = \left(\begin{array}{c} 0\\ 0 \end{array}\right) ,\,
N = S = 0.
```
Finally, notice that the corresponding transfer function reads
```math
H(s) = \dfrac{s}{ms^2+bs+k}.
```

We want to recover/identify such a structure on the basis of data only, using the LF.

### \#1.1: hand-written part

To familiarize you with LF, let us start with a very simple hand-written exercise (e.g. with $m=b=k=1$).

1. Select very simple interpolation points, e.g.
```math
\lambda = [1,2] \text{ and } \mu = -\lambda.
```
2. Evaluate the data 
```math
H(\lambda_1)=w_1,H(\lambda_2)=w_2 \text{ and }  H(\mu_1)=v_1,H(\mu_2)=v_2.
```
3. Construct the Loewner matrix $\mathbb L$.
4. Construct the shifted Loewner matrix $\mathbb M$.
5. Construct input and output data matrices (vectors) $V$ and $W$.
6. Compute the rank of $\mathbb L$.
7. Compute the eigenvalues of the matrix pencil $(\mathbb L,\mathbb M)$, being the pair $(D,V)$ solving where $D$ is a diagonal matrix with eigenvalue entries and $V$ its associated right eigenvectors (note that here, $\mathbb L$ is full column rank, thus can simplify with eigenvalues of $\mathbb L^{-1}\mathbb M$),
```math
A V = D E V.
```
8. Compute the identified rational form as
```math
R(s)= W(-s\mathbb L+\mathbb M)^{-1}V.
```
9. Conclude.

### \#1.2: numerical part

#### \#1.2.1: Define the system

Define a state-space model of the MSD and analyze it. Let us start by defining the state-space model with $m=k=b=1$. 
```Matlab
% MSD parameters
m 		= 1; % mass
k 		= 1; % spring
b 		= 1; % damper
% System properties
n   	= 2;
ny  	= 1;
nu  	= ny;
% State-space X=[x dx]
E   	= eye(n);
A   	= [0 1; -k/m -b/m];
B   	= [0; 1/m];
C   	= [0 1];
D   	= 0;
Hss 	= dss(A,B,C,D,E);
hss 	= @(s) C*((s*E-A)\B) + D;
% State-space with pH form X=ALPHA=[alpha_x alpha_v]
In  	= eye(n);
J   	= [0 1; -1 0]/m;
R   	= [0 0; 0 b]/m^2;
Q   	= [k 0; 0 m];
G   	= [0; 1/m];
P   	= [0; 0];
N   	= 0;
S   	= 0;
hph 	= @(s) ((G+P).'*Q)*((s*In-(J-R)*Q)\(G-P))+(N+S);
hph_dx 	= @(t,x,u) ((J-R)*Q)*x + (G-P)*u;
hph_y   = @(x,u) ((G+P).'*Q)*x + (N+S)*u;
```
Compute eigenvalues (singularities of $H$), zeros (zeros of $H$) and spectral zeros (zeros of $\Phi(s)=H(s)+H(-s)^\top$) of the original system.
```Matlab
eigS        = eig(A,E);
zerS        = eig([A B; C D],blkdiag(E,zeros(ny,nu)));
[SZ_v,SZ]   = lf.spectral_zeros(Hss);
```

#### \#1.2.2: First Loewner

Let us now apply the LF. One may start with the same set-up (i.e. interpolation points) as in exercise \#1.1.
```Matlab
pts     = [1 2];
la      = pts;
mu      = -pts;
k       = length(la);
q       = length(mu);
R       = ones(nu,k);
L       = ones(q,ny);
for ii = 1:k; W(1:ny,1:nu,ii) = hph(la(ii)); end
for ii = 1:q; V(1:ny,1:nu,ii) = hph(mu(ii)); end
```
Then, apply the nominal LF (tangential version) to compute an approximant:
```Matlab
[hloe,info_loe]     = lf.loewner_tng(la,mu,W,V,R,L);
```

#### \#1.2.3: Apply LF with passivity preservation

Now enforce passivity. As there is a $D$ term is equal to zero, the theorem derived in the slides (Benner et al. 2021) are no longer satisfied, use the numerical trick proposed in (Poussot-Vassal et al. 2023). And recover the pH structure:
```Matlab
opt.Ds              = 1e-2; % numerical trick to deal with non-strictly passive systems (this one is)
[hloep,info_loep]   = lf.loewner_passive(la,mu,W,V,R,L,D,opt);
[hloeph,info_loeph] = lf.passive2ph(info_loep.Hrn); % use the normalized Loewner realization
```
Analyse the outputs, such as singular values; Then, evaluate the eigenvalues, zeros, and spectral zeros...

#### \#1.2.4: Go to the time-domain simulations

If the original state-space is known, compute the equivalent Loewner projectors to recover / lift the very same states variables.
```Matlab
[Vproj,Wproj,Vproj_x0]  = lf.projectors(Hss,info_loep,'none');
```
Now, simulate in the time-domain using an exciting signal $u$ signal, e.g.
```Matlab
dt          = .01;
t           = 0:dt:20;
f           = .01; % principal harmonics
sig 		= -.1; % decay rate
u           = @(t) (sin(2*pi*f*t.^3).*exp(sig.*t))';
uu          = u(t);
x0          = zeros(n,1);
```
both full original 
```Matlab
[tt,xx]     = ode45(@(t,x) hph_dx(t,x,u(t)),t,x0);
yy          = hph_y(xx.',uu.');
```
and pH-ROM as 
```Matlab
x0r         = Vproj_x0*x0;
[tt,xr]     = ode45(@(t,x) info_loeph.dx(t,x,u(t)), tt, x0r);
yr          = info_loeph.y(xr.',uu.');
xrp         = (Vproj*xr.')';
```
Plot and compare outputs and states trajectories of the original and identified pH models, and conclude. Below a set of plots
```Matlab
mw = 10; lw = 3;
% >> Singular values
figure, hold on, grid on
plot(info_loe.sv,'-o','MarkerSize',mw,'LineWidth',lw)
plot(info_loe.sv_nu,'-o','MarkerSize',mw,'LineWidth',lw)
set(gca,'YScale','log')
xlabel('$k$','Interpreter','latex'), ylabel('Normalized singular value','Interpreter','latex')
title(['Singular values $r=' num2str(info_loe.r) '$'],'Interpreter','latex')
legend({'svd($[\bf{L},\bf{M}]$)','svd($\bf{L}$)'},'Location','best','Interpreter','latex')
% >> Eigenvalues
eigHt   = eig(info_loe.Hr);
eigHp   = eig(info_loep.Hr);
figure, hold on, grid on
plot(real(eigS),imag(eigS),'.','MarkerSize',mw,'LineWidth',lw);
plot(real(eigHt),imag(eigHt),'o','MarkerSize',mw,'LineWidth',lw);
plot(real(eigHp),imag(eigHp),'s','MarkerSize',mw,'LineWidth',lw);
title(['Eigenvalues $r=' num2str(info_loe.r) '$'],'Interpreter','latex')
xlabel('Real','Interpreter','latex'), ylabel('Imag.','Interpreter','latex')
legend({'Original' ['Loewner (\texttt{isPassive}:' num2str(isPassive(info_loe.Hr)) ')'] ['passive Loewner (\texttt{isPassive}:' num2str(isPassive(info_loep.Hr)) ')']},'location','best','Interpreter','latex')
% >> Bode magnitude
w = logspace(-2,3,1e4);
for i = 1:numel(w)
    hph_(1,1,i)   	= hph(1i*w(i));
    hloe_(1,1,i)    = hloe(1i*w(i)); 
    hloep_(1,1,i)  	= hloep(1i*w(i));
    hloeph_(1,1,i) 	= hloeph(1i*w(i));
end
figure, hold on, grid on
plot(w,20*log10(abs(squeeze(hph_(1,1,:)))),'-','LineWidth',lw), set(gca,'XScale','log')
plot(w,20*log10(abs(squeeze(hloe_(1,1,:)))),'--','LineWidth',lw), set(gca,'XScale','log')
plot(w,20*log10(abs(squeeze(hloep_(1,1,:)))),':','LineWidth',lw), set(gca,'XScale','log')
plot(w,20*log10(abs(squeeze(hloeph_(1,1,:)))),'-.','LineWidth',lw), set(gca,'XScale','log')
xlabel('pulsation [rad/s]','Interpreter','latex');
sgtitle('Bode magnitude','Interpreter','latex')
legend({'Original' ... 
       ['Loewner (\texttt{isPassive}:' num2str(isPassive(info_loe.Hr)) ')'] ...
       ['passive Loewner (\texttt{isPassive}:' num2str(isPassive(info_loep.Hr)) ')'] ... 
       'pH Loewner'},'location','best','Interpreter','latex')
% >> Time domain
figure, 
subplot(311), hold on
plot(tt,uu,'-','LineWidth',lw,'DisplayName','$u$'), grid on
legend('show','Interpreter','latex'); xlabel('$t$','Interpreter','latex'), ylabel('Input','Interpreter','latex'),
subplot(312), hold on
plot(tt,yy,'-','LineWidth',lw,'DisplayName','Original'), grid on
plot(tt,yr,'--','LineWidth',lw,'DisplayName','pH-ROM')
legend('show','Interpreter','latex'); xlabel('$t$','Interpreter','latex'), ylabel('Output','Interpreter','latex'),
subplot(313), hold on
h1=plot(tt,xx,'-','LineWidth',lw); grid on
h2=plot(tt,xr,'m:','LineWidth',lw);
h3=plot(tt,xrp,'k--','LineWidth',lw);
leg = legend([h1(1), h2(1), h3(1)], ... 
             'Original',...
             'pH-ROM', ...
             'pH-ROM (lifted)');
set(leg, 'interpreter','latex','Interpreter','latex')
xlabel('$t$','Interpreter','latex'), ylabel('Internal variables','Interpreter','latex'),
```

## Exercise \#2: First autonomous try with LF

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

Now construct the data, being the evaluation of $G$ at interpolation points $\lambda$'s and $\mu$'s (tangential directions are just set to one as the system is SISO). 
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

Examinate the outputs `htng` and `itng`, and conclude for either $G_1$ or $G_2$,
- on the approximation order $r$ (`itng.r`) and McMillan degree $\nu$ (`itng.nu`);
- on the normalized singular values of the Loewner pencil (`itng.sv`);
- on the accuracy of the eigenvalues of the obtained model;
- on frequency gain and phase response (by plotting Bode and/or Nyquist);
- on the Sylvseter equations property;
- on the informations contained in `itng` (have a deeper look at these informations).

## Exercise \#3: dig a bit more...

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

Examinate the outputs `htng_cplx` and `itng_cplx`, and conclude for either $G_1$ or $G_2$,
- on the approximation order $r$ (`itng_cplx.r`) and McMillan degree $\nu$ (`itng_cplx.nu`);
- on the normalized singular values of the Loewner pencil (`itng_cplx.sv`);
- on the accuracy of the eigenvalues of the obtained model;
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
- on the approximation order $r$ (`itng_real.r`) and McMillan degree $\nu$ (`itng_real.nu`);
- on the normalized singular values of the Loewner pencil (`itng_real.sv`);
- on the accuracy of the eigenvalues of the obtained model;
- on frequency gain and phase response (by plotting Bode and/or Nyquist);
- on the Sylvseter equations property;
- on the informations contained in `itng_real` (have a deeper look at these informations).

Notice that the matrices and values are now real.


## Exercise \#4: about reduction, passive models and pH structures

### Run the demo file `demo_ph.m`

Now we suggest to try the script `demo_ph.m`. The latter uses the very same LF but with an additional step to ensure passivity and port Hamiltonian structre. Please spend some time to test the different models proposed (or at least one).
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


## Exercise \#5: go far away and measure the flexibility of LF !

- Play with target order (`opt.target`).
- Change the original function and try replace with an irrational one such as
```math
G(s)=\frac{1}{se^{-s}+1}
```
```math
G(s)=s^3+s+1
```
- Enjoy the Loewner rule!



