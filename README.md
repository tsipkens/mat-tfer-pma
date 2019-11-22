
# MATLAB tools for PMA transfer function evaluation (mat-tfer-pma)

[![DOI](https://zenodo.org/badge/191454449.svg)](https://zenodo.org/badge/latestdoi/191454449)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)
[![Version](https://img.shields.io/badge/Version-1.4-blue.svg)]()

## 1. Code description and components

The attached MATLAB functions and scripts evaluate the transfer
function of particle mass analyzers (PMAs), including the centrifugal
particle mass analyzer (CPMA) and aerosol particle mass analyzer (APM).
This is primarily done using a novel set of expressions derived from particle
tracking methods, information for which is given in an associated paper
[(Sipkens, Olfert, and Rogak, 2019a)][ast19],
A function is also included to perform finite difference simulations of the particle number
throughout the classifier and evaluate the resultant numerical transfer function. Further information on the different methods in this distributions is given as header information in each file, such that only a brief overview is provided here.

This program contains two main components:

1. A MATLAB package, tfer_PMA, which contains the relevant functions
for evaluating the transfer function of a PMA. Such a package is
designed to be imported into other projects (e.g. it has been imported into
[mat-2d-aerosol-inversion][https://github.com/tsipkens/mat-2d-aerosol-inversion]),
where it can speed or improve the accuracy of PMA data inversion.

2. Various `main*` scripts that are used to call and analyze the transfer function
under different conditions. Of particular note is the `main` script that
evaluates the full range of available methods and produces figures similar
to those in the associate paper [(Sipkens, Olfert, and Rogak, 2019a)][ast19] and poster
[(Sipkens, Olfert and Rogak, 2019b)][eac19].

These two components are each discussed in more detail below.

## 2. The transfer function package: +tfer_pma

This package forms the core of the program and can be used
in other projects (e.g. it has been imported into
[mat-2d-aerosol-inversion][https://github.com/tsipkens/mat-2d-aerosol-inversion]).

#### 2.1 Functions to evaluate transfer functions: `tfer_*(...)`

As noted above, the core of this program is a set of
functions evaluating the transfer function for the various
cases presented in the associated work, that is
the functions featuring the names `tfer_*`. 
Further details on these functions are provided below.

###### 2.1.1 Input and output arguments

The functions share common inputs:

1. `sp` - a MATLAB structure containing the relevant setpoint parameters
(see 2.2 for the relevant function used to create this structure),

2. `m` - the masses at which the transfer function will be evaluated,

3. `d` - the mobility diameter (either as a scalar or as a vector with the
  same length as the masses at which the transfer function is to be
  evaluated),

4. `z` - the integer charge state (either as a scalar or as a vector with the
  same length as the masses at which the transfer function is to be
  evaluated),

5. `prop` - a struct that contains the properties of the particle mass analyzer
  (a sample script to generate this quantity is include as `prop_PMA.m`), and

The functions also often share common outputs:

1. `Lambda` - the transfer function and

2. `G0` - the mapping function, transforming a finial radius to the
corresponding position of the particle at the inlet (only available for
the particle tracking methods).

###### 2.1.2 Structure of filenames

Alphanumeric codes are appended to the filenames and
refer to the method or approximation used in transfer function evaluation.
Details on these alphanumeric codes are provided in the associated journal article
[(Sipkens, Olfert, and Rogak, 2019a)][ast19].
The following acts as a brief summary:

1. *FD* - Finite difference simulations to evaluate a numerical transfer function.

2. *1C* - 1st order Taylor series expansion for the particle migration velocity about the centerline radius, `rc`.

3. *1S* - 1st order Taylor series expansion for the particle migration velocity about the equilibrium radius, `rs`.

4. *2C* - 2nd order Taylor series expansion for the particle migration velocity about the centerline radius, `rc`.

5. *2S* - 2nd order Taylor series expansion for the particle migration velocity about the equilibrium radius, `rs`.

6. *W1* - Exact expression for the particle migration velocity, assuming that the inner and outer electrodes have the same rotational speed (i.e. APM conditions, `omega_hat = 1`).

7. *GE* - The exact expression for the particle migration velocity, that is a generalization of the *W1* case.

8. *ehara* - A direct implementation of the original expressions derived by [Ehara, Hagwood, and Coakley (1996)][ehara96]. This function will not be accurate for CPMA conditions (that is, it requires that `omega_hat = 1` in the `prop` structure).

9. *tri* - A triangular approximation of the transfer function using the setpoint resolution and mass of a single charged particle (multiple charging is not currently incorporated).

The codes are also occasionally modified with additional suffixes:

*pb* - Replaces the default treatment of a plug axial flow with a parabolic axial flow.

*diff* - Uses the diffusing form of the transfer function given in [(Sipkens, Olfert, and Rogak, 2019a)][ast19].

#### 2.2 Determining the setpoint: `get_setpoint(...)`

This function parses a series of name-value pairs to output a cohesive
structure fully defining the device setpoint, `sp`. This method takes
two inputs.

1. `prop` - This is the aforementioned struct that contains the properties of the
particle mass analyzer and

2. `varargin` - This is a variable length input that contains a
series of name-value pairs used to determine the setpoint. The setpoint generally
requires one to specify two of the setpoint parameters, which can include:
(a) the setpoint mass for a singly charged particle, `m_star`;
(b) the voltage, `V`;
(c) the resolution, `Rm`;
(d) the angular speed at the centerline, `omega`; and the
(e) the angular speed at the inner electrode, `omega1`;.
In general, if `m_star` is specified as one of the setpoint parameters, any one
of the other parameters can be specified. If `m_star` is not specified, the program
will expect inputs for `V` and `omega`. Other combinations are not supported. The
name-value pairs are specified as is typical of similar MATLAB functions. For example,
to specify `m_star` as 0.1 fg and `V` as 20 V, one can enter
```
sp = tfer_pma.get_setpoint(prop,'m_star',0.1e-18,'V',20);
```
If only the setpoint mass is specified as a name-value pair, the program will
use a resolution of 3 or `Rm = 3`. Accordingly,
```
sp = tfer_pma.get_setpoint(prop,'m_star',0.1e-18);
```
and
```
sp = tfer_pma.get_setpoint(prop,'m_star',0.1e-18,'Rm',3);
```
will yield identical results.

#### 2.3 Remaining functions

The remaining functions help in transfer function evaluation, with the
details provided in each file. This includes functions to convert
between particle mass and electromobility. Notably, `mp2zp.m` invokes
the mass-mobility relation to determine the mobility of particles.
There are certain assumptions implicit in this evaluation that
should be checked by the user.

## 3. Demonstration scripts: `main*`

The `main` script is included to demonstrate evaluation of the transfer function
over multiple cases. Figure 2 that is produced by this procedure will
resemble those given in the associated journal article [(Sipkens, Olfert, and Rogak, 2019a)][ast19].

Other scripts, `main_*` are intended to replicate figures in other
works and to consider multiple charging.

----------------------------------------------------------------------

#### License

This software is licensed under an MIT license (see the corresponding file
for details).

#### Contact

This program was written by Timothy A. Sipkens
([tsipkens@mail.ubc.ca](mailto:tsipkens@mail.ubc.ca)) while at the
University of British Columbia. This program was written in close
association with Steven Rogak (University of British Columbia) and
Jason Olfert (University of Alberta).

#### How to cite

This code should be cited by:

1. citing the code directly, using the DOI assigned to this code (see the archived
  versions of this code on [Zenodo](https://zenodo.org/badge/latestdoi/191454449)), and

2. citing the associated journal article describing the particle tracking methods
used in this program [(Sipkens, Olfert, and Rogak, 2019a)][ast19].

#### References

[Ehara, K., C. Hagwood, and K. J. Coakley. 1996. Novel method to classify aerosol particles according to their mass-to-charge ratio—Aerosol particle mass analyser. *J. Aerosol Sci.* 27 (2):217–34. doi: 10.1016/0021-8502(95)00562-5.][ehara96]

[Sipkens, T. A., J. S. Olfert, and S. N. Rogak. 2019a. New approaches to calculate the transfer function of particle mass analyzers. *Aerosol Sci. Technol.* doi: 10.1080/02786826.2019.1680794.][ast19]

[Sipkens, T. A., J. S. Olfert, and S. N. Rogak. 2019b. Examination of the methods available to compute the transfer function of CPMA and APM devices. Poster presented at the European Aerosol Conference. Gothenburg, Sweden, August 26.][eac19]

[ehara96]: https://doi.org/10.1016/0021-8502(95)00562-5
[ast19]: https://doi.org/10.1080/02786826.2019.1680794
[eac19]: https://www.researchgate.net/publication/336549933_Examination_of_the_methods_available_to_compute_the_transfer_function_of_CPMA_and_APM_devices
