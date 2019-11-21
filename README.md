# MATLAB tools for PMA transfer function evaluation (mat-tfer-pma)

[![DOI](https://zenodo.org/badge/191454449.svg)](https://zenodo.org/badge/latestdoi/191454449)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)

The attached MATLAB functions and scripts are intended to reproduce the 
results of the associated paper [[1][1]]. They evaluate the transfer 
function of  particle mass analyzers (PMAs), including the centrifugal 
particle mass analyzer (CPMA) and aerosol particle mass analyzer (APM). 
This is done using a novel set of expressions derived from particle 
tracking methods or using a finite difference method. Information on 
each file is given as header information in each file, and only a brief
overview is provided here. Information on the methods
is given in the associated paper [[1][1]].

## 1. Code description and components

This program contains two main components: 

1. A MATLAB package, `+tfer_PMA`, which contains the relevant functions 
for evaluating the transfer function of a PMA. Such a package is 
designed to be imported into other projects (e.g. 
https://github.com/tsipkens/mat-2d-aerosol-inversion), where it can speed or
improve the accuracy of PMA data inversion. 

2. Various `main*.m` scripts that are used to call and analyze the transfer function 
under different conditions, particularily those relevant to the associated 
paper [[1][1]] and poster [[2][2]]. 

These are each discussed in more detail below. 

## 2. The transfer function package: `+tfer_pma`

#### 2.1 Functions to evaluate transfer functions: `tfer_*.m`

As noted above, the core of this program is a set of 
functions evaluating the transfer function for the various 
cases presented in the associated work, that is
the functions featuring the names `tfer_*.m`. The file names feature case
letters, corresponding to different assumptions about the particle
migration velocity and flow conditions discussed in the associated work,
as well as `diff` and `pb` suffixes for those transfer function that
include diffusion and assume parabolic axial flow conditions, respectively.

The functions share common inputs:

1. *m_star* - the setpoint mass,

2. *m* - the masses at which the transfer function will be evaluated,

3. *d* - the mobility diameter (either as a scalar or as a vector with the
  same length as the masses at which the transfer function is to be
  evaluated),

4. *z* - the integer charge state (either as a scalar or as a vector with the
  same length as the masses at which the transfer function is to be
  evaluated),

5. *prop* - a struct that contains the properties of the particle mass analyzer
  (a sample script to generate this quantity is include as `prop_PMA.m`), and

6. *varargin* (optional) - name-value pairs to specify either the equivalent
  resolution, inner electrode angular speed, or voltage.

The functions also often share common outputs:

1. *Lambda* - the transfer function, 

2. *G0* - the mapping function, transforming a finial radius to the
corresponding position of the particle at the inlet, and

3. *sp* - a structure containing a complete description of the 
device set point used in evaluation.

Note that in these functions, there is a reference to the function
`get_setpoint.m`. This function parses the inputs *d* and *z* and then
evaluates the setpoint and related properties, including *C0*, *alpha*, and *beta*.
The key output from this function is the *sp* structure, which contains
the parameters that fully specify the PMA setpoint. 

#### 2.2 Remaining functions

The remaining functions help in transfer function evaluation, with the
details provided in each file. THis includes functions to convert
between particle mass and electromobility. Notably, `mp2zp.m` invokes
the mass-mobility relation to determine the mobility of particles. 
There are certain assumptions implicit in this evaluation that
should be checked by the user. 

## 3. Demonstration scripts: `main*.m`

The main.m script is included to demonstrate evaluation of the transfer function
over multiple cases. Figure 2 that is produced by this procedure will
resemble those given in the associated work [[1][1]]. 

Other scripts, `main_*.m` are intended to replicate figures in other 
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

#### References

1. [Sipkens et al. 2019. New approaches to calculate the transfer function of particle mass analyzers. *Aerosol Sci. Technol.* doi: 10.1080/02786826.2019.1680794.][1]
2. [Sipkens et al. 2019. Examination of the methods available to compute the transfer function of CPMA and APM devices. Poster presented at the European Aerosol Conference. Gothenburg, Sweden, August 26.][2]

[1]: https://doi.org/10.1080/02786826.2019.1680794
[2]: https://www.researchgate.net/publication/336549933_Examination_of_the_methods_available_to_compute_the_transfer_function_of_CPMA_and_APM_devices

