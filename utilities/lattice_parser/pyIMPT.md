Author: Biaobin Li

Email: biaobin@ustc.edu.cn

# Introduction

This python script could transform string-type input file into `ImpactZ.in` file. Users have experiences with `ELEGANT` will enjoy this tool. Right now, not all `IMPACT-Z` elements are added in the code. Users could go to the python scipts: `impactz_parser.py/__default_lattice()` and `impactz_parser.py/impactzin_lattice()` to add the new elements you want to use. 



# How to run it

Add `genimpactzin` into your `PATH`, for my  case, add the following line to your `.bashrc`:

```bash
export PATH=/mnt/d/githubProj/IMPACT-Z/utilities/lattice_parser:$PATH
```

Given the `lte.impz` input file:

```bash
genimpactzin lte.impz line
```

which will generate the `ImpactZ.in` file. Now you can run the `ImpactZexe` in parallel version as:

```
mpirun -np 4 ImpactZexe
```

4 processes are used as $core\_num\_T\times core\_num\_L=4$.

If you are using the single process version, `core_num_T=1,core_num_L=1` should be given. Then just type:

```
ImpactZ.exe
```

`ImpactZ.in` is automatically read.

The user is encouraged to have a look in `utilities/lattice_parser/examples`, one example is given to show how this work.



# A simple example

For the convenience of illustration, `lte.impz` refers to the python level read-in file, `ImpactZ.in` refers to `ImpactZexe` read-in file. You can rename `lte.impz`  any other names you like.

The `lte.impz` file mainly consists of three sections, `control, beam and lattice` sections. The detailed mapping relationships between`ImpactZ.in` and `lte.impz` are listed in the following section. Here we give a simple example for the  usage of `lte.impz` :

```python
!control section
!===============
&control

core_num_T = 2;
core_num_L = 2;
meshx = 32;
meshy = 32;
meshz = 64;
kinetic_energy = 300e6;
freq_scale = 1.3e9;

&end

!beam section
!==============
&beam
mass = 0.511001e6;
charge = -1.0;

distribution_type = 2;
Np = 5000;
total_charge = 1e-9;

emit_nx=0.176e-6, beta_x=12.73, alpha_x=-0.85;
emit_ny=0.176e-6, beta_y=12.73, alpha_y=-0.85;

sigz=1e-3, sigE=5e3;

&end

!lattice section
!=====================
&lattice

!rpn expression is supported,
!only a few mathematical operator are added, please see 
!lattice_parser.py/rpn_cal() for more details.
!------------------------------------------------------
% 0.2 sto LB1
% -4.410 pi * 180 / sto AB1  ! Bend angle

BCX11: BEND,L= LB1,ANGLE=AB1,       E2=AB1,       steps=1, pipe_radius=2.1640E-02, fint=0.3893
BCX12: BEND,L= LB1,ANGLE= "0 AB1 -",E1= "0 AB1 -",steps=1, pipe_radius=2.1640E-02, fint=0.3893
BCX13: BEND,L= LB1,ANGLE= "0 AB1 -",E2= "0 AB1 -",steps=1, pipe_radius=2.1640E-02, fint=0.3893
BCX14: BEND,L= LB1,ANGLE=AB1,       E1=AB1,       steps=1, pipe_radius=2.1640E-02, fint=0.3893

D1  : DRIF, L=5.0
Dm  : DRIF, L=0.5

BC1 : LINE=(BCX11,D1,BCX12,Dm,BCX13,D1,BCX14)

W0:   watch, filename_ID=1000
W1:   watch, filename_ID=1001

Line : LINE=(W0,BC1,W1)

&end
```



# Control and beam section

The mapping relationship between `ImpactZ.in` and `lte.impz` in `control` sections are listed as following:

```bash
line1: 
core_num_T core_num_L

line2:
6 Np integrator error output_ratio

line3：
meshx meshy meshz flagbc x_pipe_width y_pipe_width period_len

line4:
distribution_type restart sub_cycle 1

line5:
Np

line6:
current  #Q*f_scale

line7:
# value defined automatically by charge/mass in the code; 
# the definition of charge and mass, see beam section.

line8-line10:
# defined by beam section
alpha_x beta_x emit_x mismatchx mismatchpx offsetX     offsetPx
alpha_y beta_y emit_y mismatchy mismatchpy offsetY     offsetPy
alpha_z beta_z emit_z mismatchz mismatchE  offsetPhase offsetEnergy

line11:
current kinetic_energy mass charge freq_scale ini_phase 

```



## Control parameters

All control parameters in `lte.impz` are listed:



| Parameter Name | Units | Type   | Default | Description                                                  |
| -------------- | ----- | ------ | ------- | ------------------------------------------------------------ |
| core_num_T     |       | int    | 1       | processor number for the transverse direction.               |
| core_num_L     |       | int    | 1       | processor number for the longitudinal direction.             |
| dt             | s     | double | 1e-12   | time step size                                               |
| max_step       |       | int    | 1e6     | maximum number of time steps                                 |
| nbunch         |       | int    | 1       | see manual for details.                                      |
| Dim            |       | int    | 6       | random seed integer                                          |
|                |       |        |         |                                                              |
|                |       |        |         |                                                              |
| error          |       | int    | 0       | Error study?                                                 |
| diag           |       | int    | 1       | see manual for details.                                      |
| image_sc       |       | int    | 1       | Image charge flag.                                           |
| image_stop_pos | m     | double | 0.02    | position image charge forces are neglected.                  |
| meshx          |       | int    | 32      | number of mesh points in x direction.                        |
| meshy          |       | int    | 32      | number of mesh points in y direction.                        |
| meshz          |       | int    | 32      | number of mesh points in z direction.                        |
|                |       |        |         |                                                              |
| Xrad           | m     | double | 0.015   | size of computational domain. Transverse size.               |
| Yrad           | m     | double | 0.015   | size of computational domain. Transverse size.               |
| PerdLen        | m     | double | 10.0    | PerdLen should be greater than the beam line lattice length. |
|                |       |        |         |                                                              |
| Restart        |       | int    | 0       |                                                              |
| Nemission      |       | int    | 400     | the number of numerical emission steps.                      |
| Temission      | s     | double | 1e-9    | Laser pulse emission time.                                   |
| kinetic_energy | eV    | double | 0       | The real kinetic energy of the beam. NOT Bkenergy in line9.  |
| freq_rf_scale  | Hz    | double | 2856e6  | scale frequency $f_{scal}$,  $Scxl=c/(2\pi f_{scal}) $.      |
| ini_t          | s     | double | 0.0     | initial reference time.                                      |



## Beam parameters

All beam section parameters in `lte.impz` are listed:

| Parameter Name    | Units | Type   | Default | Description                                                |
| ----------------- | ----- | ------ | ------- | ---------------------------------------------------------- |
| mass              | eV    | double | 0.511e6 | mass of the particle.                                      |
| charge            |       | double | -1.0    | -1 for electron.                                           |
| distribution_type |       | int    | 2       | 6D gaussian distribution. See more options in Ji’s manual. |
| Np                |       | int    | 1e3     | particle number.                                           |
| total_charge      | C     | double | 1e-9    | charge of the beam.                                        |
|                   |       |        |         |                                                            |
| emit_x            | m rad | double | 0.0     | emitance.                                                  |
| emit_nx           | m rad | double | 0.0     | normalized emittance.                                      |
| beta_x            | m     | double | 1.0     | twiss para.                                                |
| alpha_x           |       | double | 0.0     | twiss para.                                                |
| sigx              | m     | double | 0.0     | rms bunch size.                                            |
| sigpx             |       | double | 0.0     | rms value of $\gamma\beta_x/\gamma_0\beta_0$               |
|                   |       |        |         |                                                            |
| emit_y            | m rad | double | 0.0     | emittance.                                                 |
| emit_ny           | m rad | double | 0.0     | normalized emittance.                                      |
| beta_y            | m     | double | 1.0     | twiss para.                                                |
| alpha_y           |       | double | 0.0     | twiss para.                                                |
| sigy              | m     | double | 0.0     | rms bunch size.                                            |
| sigpy             |       | double | 0.0     | rms value of $\gamma\beta_y/\gamma_0\beta_0$               |
|                   |       |        |         |                                                            |
| emit_z            | m rad | double | 0.0     | emittance.                                                 |
| emit_nz           | m rad | double | 0.0     | normalized emittance.                                      |
| beta_z            | m     | double | 1.0     | twiss para.                                                |
| alpha_z           |       | double | 0.0     | twiss para.                                                |
| sigz              | m     | double | 0.0     | rms bunch length.                                          |
| sigpz             | eV    | double | 0.0     | rms value of $\gamma\beta_z/\gamma_0\beta_0$               |
|                   |       |        |         |                                                            |

Users could either use twiss parameters to define initial beam distribution, or use rms values. For $\sigma_{ij}\neq0$ cases, please use twiss-para. 



In the definition: $z=ct$.



# Lattice section

Right now, only a few frequently used elements in `ImpactZ.in` are added into the python parser.



## Elements

### DRIFT

0 element.

| Parameter Name | Units | Type   | Default | Description     |
| -------------- | ----- | ------ | ------- | --------------- |
| zedge          | m     | double | 0.0     | global position |
| L              | m     | double | 0.0     | length of drift |



### QUAD

1 element.

| Parameter Name | Units | Type       | Default | Description                                                  |
| -------------- | ----- | ---------- | ------- | ------------------------------------------------------------ |
| zedge          | m     | double     | 0.0     | global position                                              |
| L              | m     | double     | 0.0     | length                                                                             |
| grad           | T/m   | double     | 0.0     | quadrupole strength, $gradient=\frac{\partial B_y}{\partial x}$ |
| fileid        |       | int/double | 1.0     | use minus x to refer rfdatax; when fileid=1.0, use Enge function for fringe field, fileid is changed to quad length in python scripts. |
|                |       |            |         |                                                              |
| Dx             | m     | double     | 0.0     | x misalignment error                                         |
| Dy             | m     | double     | 0.0     | y misalignment error                                         |
| rotate_x       | rad   | double     | 0.0     | rotation error in x direction                                |
| rotate_y       | rad   | double     | 0.0     | rotation error in y direction                                |
| ratate_z       | rad   | double     | 0.0     | rotation error in y direction                                |
| freq           | Hz    | double     | 0.0     | rf quadrupole frequency                                      |
| phase          | deg   | double     | 0.0     | rf quadrupole phase                                          |



### SOL

3 element.

| Parameter Name | Units | Type   | Default | Description     |
| -------------- | ----- | ------ | ------- | --------------- |
| zedge          | m     | double | 0.0     | global position |
| L              | m     | double | 0.0     | length          |
| fileid        |       | int    | None    | file ID         |



### SOLRF

105 element.

| Parameter Name | Units | Type   | Default | Description                   |
| -------------- | ----- | ------ | ------- | ----------------------------- |
| zedge          | m     | double | 0.0     | global position.              |
| L | m | double | 0.0 | the longitudinal length of the element (Blength). |
| Emax       | V/m | double | 1.0     | the absolute maximum values of on-axis Ez field. |
| freq           | Hz    | double | 2856e6  | frequency                     |
| phase          | deg   | double | 0.0     | RF design phase               |
| Dx             | m     | double | 0.0     | x misalignment error          |
| Dy             | m     | double | 0.0     | y misalignment error          |
| rotate_x       | rad   | double | 0.0     | rotation error in x direction |
| rotate_y       | rad   | double | 0.0     | rotation error in y direction |
| ratate_z       | rad   | double | 0.0     | rotation error in y direction |
| scaleB         |       | double | 0.0     | scale of solenoid B field.    |
|  | |  |  |  |
| fileid | | int | None | file ID |
| z1 | m | double | None | rfdatax second line. Distance before the zedge. `None`, no update. |
| z2 | m | double | None | rfdatax third line. Distance after the zedge. `None`, no update. |
| L_fourier_exp | m | double | None | rfdatax fourth line. The length of the reconstructed field using the fourier coefficients given in rfdatax. (z1,z2,L_fourier_exp) will be used to update rfdatax line2-line4. `None`, no update. |







### WATCH

-2 element.

| Parameter Name | Units | Type   | Default | Description      |
| -------------- | ----- | ------ | ------- | ---------------- |
| zedge          | m     | double | 0.0     | global position  |
|                |       |        |         |                  |
| filename_id    |       | int    | 80      | fort.80          |
| sample_freq    |       | int    | 1       | sample out freq. |





### TWS

Traveling wave structure, without entrance and exit coupler included.

| Parameter Name | Units | Type   | Default  | Description                                                  |
| -------------- | ----- | ------ | -------- | ------------------------------------------------------------ |
| zedge          | m     | double | 0.0      | global position.                                             |
| Lperd          | m     | double | 0.105    | one period length of the field, for $2\pi/3$ cavity, `Lperd` is the length of 3 cells. |
| Nperd          |       | double | 1.0      | how many periods to repeat. `Blength` and 4th line in `rfdatax` will be updated with $Lperd\times Nperd$. |
| phaseshift     | rad   | double | $2\pi/3$ | phase shift of one cell.                                     |
| Emax           | V/m   | double | 0.0      | the absolute maximum values of on-axis Ez field.             |
| freq           | Hz    | double | 2856e6   | RF frequency                                                 |
| phase          | deg   | double | 0.0      | RF design phase, it should be same with the entrance coupler. |
| fileid_1       |       | int    | None     | file ID for the standing wave.                               |
| fileid_2       |       | int    | None     | file ID for the shifted standing wave.                       |
| Dx             | m     | double | 0.0      | x misalignment error                                         |
| Dy             | m     | double | 0.0      | y misalignment error                                         |
| rotate_x       | rad   | double | 0.0      | rotation error in x direction                                |
| rotate_y       | rad   | double | 0.0      | rotation error in y direction                                |
| ratate_z       | rad   | double | 0.0      | rotation error in y direction                                |

The traveling wave structure is modeled by two standing wave, one should use `RFcoeflcls` to get the fourier coefficients of the standing wave and the shifted standing wave, i.e. `rfdatafieldid_1` and `rfdatafileid_2`. Only profile information are given for the fourier coefficients.



The `z1,z2` and `length` in rfdatax will be updated according to `Lcell` and `Ncell`.
