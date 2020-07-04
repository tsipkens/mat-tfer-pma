
# Matlab tools for PMA transfer function evaluation (mat-tfer-pma)

[![DOI](https://zenodo.org/badge/191454449.svg)](https://zenodo.org/badge/latestdoi/191454449)
[![MIT license](https://img.shields.io/badge/License-MIT-blue.svg)](https://lbesson.mit-license.org/)
[![Version](https://img.shields.io/badge/Version-2.0.1-blue.svg)]()

**Note**: A python version of some of these transfer functions is available at https://github.com/tsipkens/py-tfer-pma. 

The attached Matlab functions evaluate the transfer function of particle mass analyzers (PMAs), including the centrifugal particle mass analyzer (CPMA) and aerosol particle mass analyzer (APM). This is primarily done using a novel set of expressions derived from particle tracking methods, information for which is given in an associated paper [(Sipkens, Olfert, and Rogak, 2020a)][ast20]. A function is also included to perform finite difference simulations of the particle number throughout the classifier and evaluate the resultant numerical transfer function. Further information on the different methods in this program is given as header information in each file, such that only a brief overview is provided here.

This program contains two main components:

1. The upper directory contains the relevant functions for evaluating the transfer function of a PMA. A good demonstration of how these functions can be imported into other projects is [mat-2d-aerosol-inversion](https://github.com/tsipkens/mat-2d-aerosol-inversion),
where it is used to speed or improve the accuracy of PMA data inversion.

2. The `test/`  folder contains a series of `main*` scripts in, which are used to call and analyze the transfer function under different conditions. Of particular note is the `main` script that evaluates the full range of available methods and produces figures similar to those in the associate paper [(Sipkens, Olfert, and Rogak, 2020a)][ast20] and poster [(Sipkens, Olfert and Rogak, 2019b)][eac19].

These two components are each discussed in more detail below.

## 1. Upper directory and evaluating the transfer function

The functions in the upper directory general fall into three categories.

### 1.1 Methods to evaluate transfer functions: tfer_*(...)

The  `tfer_*` functions form the core of the program and evaluate the transfer function for the various cases presented in the associated work [(Sipkens, Olfert, and Rogak, 2020a)][ast20].

Alphanumeric codes are appended to the filenames and refer to the method or approximation used in transfer function evaluation. Details on these alphanumeric codes are provided in the associated journal article [(Sipkens, Olfert, and Rogak, 2019a)][ast20]. The following acts only as a brief summary:

| Code  | Description                                                  |
| :---- | :----------------------------------------------------------- |
| FD    | Finite difference simulations to evaluate a numerical transfer function. |
| 1C    | 1st order Taylor series expansion for the particle migration velocity about the centerline radius, `rc`. |
| 1S    | 1st order Taylor series expansion for the particle migration velocity about the equilibrium radius, `rs`. |
| 2C    | 2nd order Taylor series expansion for the particle migration velocity about the centerline radius, `rc`. |
| 2S    | 2nd order Taylor series expansion for the particle migration velocity about the equilibrium radius, `rs`. |
| W1    | Exact expression for the particle migration velocity, assuming that the inner and outer electrodes have the same rotational speed (i.e. APM conditions, `omega_hat = 1`). |
| GE    | The exact expression for the particle migration velocity, that is a generalization of the *W1* case. |
| ehara | A direct implementation of the original expressions derived by [Ehara, Hagwood, and Coakley (1996)][ehara96]. This function will not be accurate for CPMA conditions (that is, it requires that `omega_hat = 1` in the `prop` structure). |
| tri   | A triangular approximation of the transfer function using the setpoint resolution and mass of a single charged particle (multiple charging is not currently incorporated). |

The codes are also occasionally modified with additional suffixes with the following meanings:

| Code | Description                                                  |
| :--- | :----------------------------------------------------------- |
| pb   | Replaces the default treatment of a plug axial flow with a parabolic axial flow. |
| diff | Uses the diffusing form of the transfer function given in [(Sipkens, Olfert, and Rogak, 2019a)][ast20]. |

##### Input and output arguments

The  `tfer_*` functions share common inputs:

1. `sp` - A Matlab structure containing the relevant setpoint parameters (see [1.2](1.2-determining-the-setpoint:-get_setpoint(...)) for the relevant function used to create this structure).

2. `m` - The masses at which the transfer function will be evaluated.

3. `d` - The mobility diameters (either as a scalar or as a vector with the same length as `m`) at which the transfer function is to be evaluated).

4. `z` - The integer charge state (either as a scalar or as a vector with the same length as `m`) at which the transfer function is to be evaluated.

5. `prop` - A Matlab structure that contains the properties of the particle mass analyzer. This includes fields that describe the inner radius, `r1`; outer radius, `r2`; length, `L`; and flows, `Q*`. A sample script to generate this quantity is included as `prop_pma(...)` in the upper directory.

The functions also often share common outputs:

1. `Lambda` - The transfer function.

2. `G0` - The mapping function, transforming a finial radius to the corresponding position of the particle at the inlet (only available for the particle tracking methods).

### 1.2 Determining the setpoint: get_setpoint(...)

This function parses a series of name-value pairs to output a cohesive structure fully defining the device setpoint, output in the form of `sp` or a *setpoint structure*. A sample set of values for this structure is found below (the precise values will depend on the PMA properties and mass-mobility relation values set in the device):

| Variable | *m*<sup>\*</sup> [kg] | *V* [V] | *R*<sub>m</sub> | **ω** [rad/s] | **ω**<sub>1 </sub> [rad/s] | **ω**<sub>2</sub> [rad/s] | **α** [s<sup>-1</sup>] | **β** [m<sup>2</sup>⋅s<sup>-1</sup>] | *m*<sub>max</sub> [kg] |
| ------ | :------------------------: | :-------: | :--------------------: | :--------------: | :--------------------------: | :--------------------------: | :-------------: | :------------: | :-------------------------: |
| Field name | `m_star` | `V` | `Rm` | `omega1` | `omega1` | `omega2` | `alpha` | `beta` | `m_max` |
| 1      |    1.0×10<sup>-20</sup>    |   24.4    |           10           |     2,583.9      |           2,583.5            |           2,505.2            |       176       |      8.66      |    1.1×10<sup>-20</sup>     |
| 2      |    1.0×10<sup>-19</sup>    |    104    |           10           |     1,659.5      |           1,684.3            |           1,633.3            |       115       |      5.65      |    1.1×10<sup>-19</sup>     |
| 3      |    1.0×10<sup>-18</sup>    |    398    |           10           |     1,026.5      |           1,042.4            |           1,010.8            |      71.0       |      3.50      |    1.1×10<sup>-18</sup>     |
| 4      |    1.0×10<sup>-17</sup>    |   1,280   |           10           |      582.69      |            591.77            |            573.83            |      40.3       |      1.99      |    1.1×10<sup>-17</sup>     |
| ...    |                            |           |                        |                  |                              |                              |                 |                |                             |

The function itself takes two inputs:

1. `prop` - This is the aforementioned struct that contains the properties of the particle mass analyzer, such as its dimensions and the flow rate, and

2. `varargin` - This is a variable length input that contains a series of name-value pairs used to determine the setpoint.

For the latter quantity, the setpoint generally requires the practitioner to specify two of the setpoint parameters, which here can include:

1. the setpoint mass for a singly charged particle, `m_star`;
2. the voltage, `V`;
3. the resolution, `Rm`;
4. the angular speed at the centerline, `omega`; and
5. the angular speed at the inner electrode, `omega1`;.

For this program, if `m_star` is specified as one of the setpoint parameters, any one of the other parameters can be specified. If `m_star` is not specified, the program will expect inputs for `V` and `omega`. Other combinations are not currently supported. The name-value pairs are specified similar to other Matlab functions. For example, to specify `m_star` as 0.1 fg and `V` as 20 V, one can enter

```Matlab
sp = get_setpoint(prop,'m_star',0.1e-18,'V',20);
```
If only the setpoint mass is specified as a name-value pair, the program will
use a resolution of 3 or `Rm = 3`. Accordingly,
```Matlab
sp = get_setpoint(prop,'m_star',0.1e-18);
```
and
```Matlab
sp = get_setpoint(prop,'m_star',0.1e-18,'Rm',3);
```
will yield identical results. The function can also handle vector inputs and will output a structured array with one entry per setpoint. For example, a vector of mass setpoints and a resolution of *R*<sub>m</sub> = 10 can be specified using:

```Matlab
m_star = 1e-18.*logspace(log10(0.1),log10(100),25); % mass setpoints
sp = get_setpoint(prop,...
	'm_star',m_star,'Rm',10); % PMA setpoints for Rm = 10
```

Note that the input to the function must either be (a) two vectors of the  same length or (b) a scalar and a vector (as in the example above).

#### 1.3 Remaining functions

The remaining functions help in transfer function evaluation, with the details provided in each file. This includes functions to convert between particle mass and electromobility.

Notably, `mp2zp.m` invokes the mass-mobility relation to determine the mobility of particles. There are certain assumptions implicit in this evaluation that should be checked by the user. The mass-mobility exponent, `Dm`, and effective density, `rho0`, used when invoking the mass-mobility relation are specified in the `prop` structure. The mobility diameter is then estimated as:

```Matlab
d = (m./prop.rho0).^(1/prop.Dm);
```

The default values can be found in that function (`prop.Dm = 3;` and `prop.rho0 = 1000*pi/6;` at the time of writing).  

## 2. Test scripts: main*

The scripts in the `test/` folder are included to demonstrate evaluation of the transfer function over multiple cases. These scripts first initialize the PMA setpoint, `sp`; properties, `prop`; and vectors at which the transfer function will be evaluated, `m`, `d` and `z`.  The scripts proceed by evaluating the transfer function and plotting the results.

In order to use the scripts, one must first add the `test/` folder to the Matlab path. When in the `mat-tfer-pma` directory, enter

```Matlab
addpath test;
```
on the MATLAB command line. The scripts can then be called by (i) entering their name on the command line or (ii) expanding the folder in the file explorer within Matlab, opening the file, and pressing "Run" (i.e., the green play button) in the Matlab ribbon.  

The `main` script, without any other text appended, is included to replicate the results of [(Sipkens, Olfert, and Rogak, 2020a)][ast20], where Figure 2 that is produced by this procedure will resemble the figures in that article.

Other scripts, `main_*` are intended to replicate figures in other works and to consider multiple charging.

## 3. Using this repository in other projects

This program is designed to be used in other projects. The [mat-2d-aerosol-inversion](https://github.com/tsipkens/mat-2d-aerosol-inversion) project, for example, uses this project as a submodule and greatly speeds the generation of kernel required for 2D inversions, as per [Sipkens, Olfert, and Rogak, (2020b)][jas20]).

----------------------------------------------------------------------

#### License

This software is licensed under an MIT license (see the corresponding file or details).

#### Contact

This program was written by Timothy A. Sipkens ([tsipkens@mail.ubc.ca](mailto:tsipkens@mail.ubc.ca)) while at the University of British Columbia. This program was written in close association with Steven Rogak (University of British Columbia) and Jason Olfert (University of Alberta).

#### How to cite

This code should be cited by:

1. citing the associated journal article describing the particle tracking methods used in this program [(Sipkens, Olfert, and Rogak, 2019a)][ast20], and

2. citing the code directly (either using the DOI assigned to the version of code used - see the archived versions of this code on [Zenodo](https://zenodo.org/badge/latestdoi/191454449) - or, less formally, making reference to the GitHub repository at https://github.com/tsipkens/mat-tfer-pma).

#### References

[Ehara, K., C. Hagwood, and K. J. Coakley. 1996. Novel method to classify aerosol particles according to their mass-to-charge ratio—Aerosol particle mass analyser. *J. Aerosol Sci.* 27:2, 217–34. DOI: 10.1016/0021-8502(95)00562-5.][ehara96]

[Sipkens, T. A., J. S. Olfert, and S. N. Rogak. 2019. Examination of the methods available to compute the transfer function of CPMA and APM devices. Poster presented at the European Aerosol Conference. Gothenburg, Sweden, August 26.][eac19]

[Sipkens, T. A., J. S. Olfert, and S. N. Rogak. 2020a. New approaches to calculate the transfer function of particle mass analyzers. *Aerosol Sci. Technol.* 54:1, 111-127. DOI: 10.1080/02786826.2019.1680794.][ast20]

[Sipkens, T. A., Olfert, J. S., & Rogak, S. N. 2020b. Inversion methods to determine two-dimensional aerosol mass-mobility distributions: A critical comparison of established methods. *J. Aerosol Sci.* 104, 105484. DOI: 10.1016/j.jaerosci.2019.105484][jas20]

[ehara96]: https://doi.org/10.1016/0021-8502(95)00562-5
[ast20]: https://doi.org/10.1080/02786826.2019.1680794
[eac19]: https://www.researchgate.net/publication/336549933_Examination_of_the_methods_available_to_compute_the_transfer_function_of_CPMA_and_APM_devices
[jas20]: https://doi.org/10.1016/j.jaerosci.2019.105484
