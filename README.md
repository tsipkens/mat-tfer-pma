## UBC-tfer-PMA

The attached functions and script are intended to reproduce the results of
the associated paper (submitted). They evaluate the transfer function of
the centrifugal particle mass analyzer (CPMA) and aerosol particle mass
analyzer (APM). This is done using a novel set of expressions derived from
particle tracking methods and using a finite difference method. Information
on each file is given as header information in each file, and only a brief
overview is provided here.

----------------------------------------------------------------------

### Code description and components

#### Functions to evaluate transfer functions (`tfer_*.m`)

The core of this program is a set of functions evaluating the transfer
function for the various cases presented in the associated work, that is
the functions featuring the names `tfer_*.m`. The file names feature case
letters, corresponding to different assumptions about the particle
migration velocity and flow conditions discussed in the associated work,
as well as `diff` and `pb` suffixes for those transfer function that
include diffusion and assume parabolic axial flow conditions, respectively.

The functions share common inputs:

1. m_star - the setpoint mass,

2. m - the masses at which the transfer function will be evaluated,

3. d - the mobility diameter (either as a scalar or as a vector with the
  same length as the masses at which the transfer function is to be
  evaluated),

4. z - the integer charge state (either as a scalar or as a vector with the
  same length as the masses at which the transfer function is to be
  evaluated),

5. prop - a struct that contains the properties of the particle mass analyzer
  (a sample script to generate this quantity is include as `prop_CPMA.m`), and

6. varargin (optional) - name-value pairs to specify either the equivalent
  resolution, inner electrode angular speed, or voltage.

Note that in these functions, there is a reference to the script
`get_setpoint.m`. This script parses the d and z inputs to each of the
`tfer_*/m` functions and then evaluates the setpoint and related properties,
including C0, alpha, and beta. 

#### Demonstration script (`main.m`)

This script is included to demonstrate evaluation of the transfer function
over multiple cases.

Figure 2 that is produced by this procedure will resemble those given in
the associated work

#### Remaining functions

The remaining functions help in transfer function evaluation, with the
details provided in each file.

----------------------------------------------------------------------

#### License

This code is distributed under an MIT license
(see corresponding LICENSE file).

----------------------------------------------------------------------

#### Contact

This program was written by Timothy A. Sipkens
([tsipkens@mail.ubc.ca](mailto:tsipkens@mail.ubc.ca)) while at the
University of British Columbia
