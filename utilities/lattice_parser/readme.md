Author: Biaobin Li

Email: biaobin@ustc.edu.cn



## Control parameters

All control parameters in `lte.impt` are listed:



| Parameter Name | Units | Type   | Default | Description                                                  |
| -------------- | ----- | ------ | ------- | ------------------------------------------------------------ |
| core_num_T     |       | int    | 1       | processor number for the transverse direction.               |
| core_num_L     |       | int    | 1       | processor number for the longitudinal direction.             |
| dt             | s     | double | 1e-12   | time step size                                               |
| max_step       |       | int    | 1e6     | maximum number of time steps                                 |
| nbunch         |       | int    | 1       | see manual for details.                                      |
| Dim            |       | int    | 6       | random seed integer                                          |
| error          |       | int    | 0       | Error study?                                                 |
| diag           |       | int    | 1       | see manual for details.                                      |
| image_sc       |       | int    | 0       | 0/1, Image charge OFF or ON.                                 |
| image_stop_pos | m     | double | 0.02    | position image charge forces are neglected.                  |
| meshx          |       | int    | 32      | number of mesh points in x direction.                        |
| meshy          |       | int    | 32      | number of mesh points in y direction.                        |
| meshz          |       | int    | 32      | number of mesh points in z direction.                        |
| Xrad           | m     | double | 0.015   | size of computational domain. Transverse size.               |
| Yrad           | m     | double | 0.015   | size of computational domain. Transverse size.               |
| PerdLen        | m     | double | 10.0    | PerdLen should be greater than the beam line lattice length. |
| Restart        |       | int    | 0       |                                                              |
| Nemission      |       | int    | -1      | the number of numerical emission steps.  `Nemission=-1` no cathode model.`Nemission=400` means it will take 400 steps for emission process, and the emission time step is `Temission/Nemission`. |
| Temission      | s     | double | 0.0     | Laser pulse emission time. `Temission=1e-9` means laser pulse is 1ns. |
| kinetic_energy | eV    | double | 0       | The real kinetic energy of the beam.                         |
| Bkenergy       | eV    | double | None    | If None, `kinetic_energy` values will be set in the python level. Particles behind the cathod will use this beta for emission. ==Question==: Under what situation would this value differs with `kinetic_energy`? |
| freq_rf_scale  | Hz    | double | 2856e6  | scale frequency $f_{scal}$,  $Scxl=c/(2\pi f_{scal}) $.      |
| ini_t          | s     | double | 0.0     | initial reference time.                                      |



`Nemission` is the total steps for emission process, `Temission` is the emission time. 

```
dt=1e-12
max_step=1e3

Nemission=100
Temission=1e-9
```

Then the emission process of 1ns would only takes 100 steps, and another 900 steps for `dt=1e-12`. When `Nemission=-1`, no cathode model, The particles are assumed to start in a vacuum. The emission process is also `dt=1e-12s`, and emission steps are 1000. 



## Beam parameters

All beam section parameters in `lte.impt` are listed:

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
| dx                | m     | double | 0.0     | offset for x                                               |
| emit_y            | m rad | double | 0.0     | emittance.                                                 |
| emit_ny           | m rad | double | 0.0     | normalized emittance.                                      |
| beta_y            | m     | double | 1.0     | twiss para.                                                |
| alpha_y           |       | double | 0.0     | twiss para.                                                |
| sigy              | m     | double | 0.0     | rms bunch size.                                            |
| sigpy             |       | double | 0.0     | rms value of $\gamma\beta_y/\gamma_0\beta_0$               |
| dy                | m     | double | 0.0     | offset for y                                               |
| emit_z            | m rad | double | 0.0     | emittance.                                                 |
| emit_nz           | m rad | double | 0.0     | normalized emittance.                                      |
| beta_z            | m     | double | 1.0     | twiss para.                                                |
| alpha_z           |       | double | 0.0     | twiss para.                                                |
| sigz              | m     | double | 0.0     | rms bunch length.                                          |
| sigpz             | eV    | double | 0.0     | rms value of $\gamma\beta_z/\gamma_0\beta_0$               |
| dz                | m     | double | 0.0     | offset for z                                               |

Users could either use twiss parameters to define initial beam distribution, or use rms values. For $\sigma_{ij}\neq0$ cases, please use twiss-para. 

In the definition: $z=ct$.

For ijk, like `distribution_type=112`, the `sigx,sigy` actually is beam radius `r`, and `sigz` is `Lbunch` full length, which in transverse direction is circle uniform, in longitudinal direction is flat-top, and z is in $[-Lbunch, 0]$ range. `zscale` is automatically given `1e-9` value in the python code as the manual said.



### Distribution

#### 112, cylinder uniform

112 could be used to generate cylinder uniform distribution.

$\sigma_{i,j,k}$ are actually refers to:
$$
sigx=sigy \\
r\equiv sigx \\
L_{bunch}\equiv sigz
$$




#### 6, cylinder uniform

New added:

The following transformation is done in the Fortran source code:
$$
\sigma_x=\sigma_y \\
r=2\times \sigma_x \\
L_{bunch}=2\sqrt{3}\sigma_z
$$
So actually used definition of the input paras. are:
$$
sigx=sigy \\
r\equiv sigx \\
L_{bunch}\equiv sigz
$$




# Lattice section

Right now, only a few frequently used elements in `ImpactT.in` are added into the python parser.



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

For `fileid=3`, the B field file `1T3.T7`, unit is [cm] and [gauss]. It is

`1T3.T7out` would be generated, five columns data:

```bash
s[m], Br(r=0) [gauss], Br(r=0+dr) [gauss], Bz(r=0) [gauss], Bz(r=dr) [gauss]
```







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
| fileid | | int | None | file ID |
| z1 | m | double | None | rfdatax second line. Distance before the zedge. `None`, no update. |
| z2 | m | double | None | rfdatax third line. Distance after the zedge. `None`, no update. |
| L_fourier_exp | m | double | None | rfdatax fourth line. The length of the reconstructed field using the fourier coefficients given in rfdatax. (z1,z2,L_fourier_exp) will be used to update rfdatax line2-line4. `None`, no update. |

The traveling wave structure is modeled by two standing wave, one should use `RFcoeflcls` to get the fourier coefficients of the standing wave and the shifted standing wave, i.e. `rfdatax`. Only profile information are given for the fourier coefficients. 

The `z1,z2` and `length` in rfdatax will be updated according to `Lcell` and `Ncell`.



### TWS

Traveling wave structure. 

| Parameter Name | Units | Type   | Default  | Description                                                  |
| -------------- | ----- | ------ | -------- | ------------------------------------------------------------ |
| zedge          | m     | double | 0.0      | global position.                                             |
| L              | m     | double | 0.0      | length of the whole TWS structure linac with couplers excluded. L should be $n\times Lcav$ |
| cav_mode       | rad   | double | $2\pi/3$ | to get $Amp=Emax/sin(\beta_0d)$                              |
| Lcoup          | m     | double | 0.052464 | entrance and exit coupler field length. Default is SLAC-S-band values. To set the global position, rfdatax will also be updated. |
| Lcav           | m     | double | 0.104926 | one period length of the trwave. For $2\pi/3$  mode S-band, `Lcav=0.104926 m`. rfdatax will be updated using this value. |
| Emax           | V/m   | double | 0.0      | the absolute maximum values of on-axis Ez field.             |
| freq           | Hz    | double | 2856e6   | RF frequency                                                 |
| phase          | deg   | double | 0.0      | RF design phase.                                             |
| fileid_1       |       | int    | None     | file ID for the entrance coupler.                            |
| fileid_2       |       | int    | None     | Auto update to `fileid_1+1` in python code.                  |
| fileid_3       |       | int    | None     | Auto update to `fileid_1+2 in python code.                   |
| fileid_4       |       | int    | None     | Auto update to `fileid_1+3 in python code. file ID for the exit coupler. |
| Dx             | m     | double | 0.0      | x misalignment error                                         |
| Dy             | m     | double | 0.0      | y misalignment error                                         |
| rotate_x       | rad   | double | 0.0      | rotation error in x direction                                |
| rotate_y       | rad   | double | 0.0      | rotation error in y direction                                |
| ratate_z       | rad   | double | 0.0      | rotation error in y direction                                |

usage:

```
tws: tws, zedge=0.0, L=3*0.104926,Emax=25.5e6,freq=2856e6,phase=0.0,fileid_1=4, Lcoup=5.2464e-2,Lcav=0.104926
```

rfdata4, rfdata5, rfdata6, rfdata7 should be given.



### EMFLDCYL

112 element. Read in discrete EM field data in cylinder coordinates. I changed the source code,  the data formats are:

```
open(33,file=’1T1.T7’)
write(33,*)rmin,rmax,nmr-1
write(33,*)zmin,zmax,nmz-1
do j=1,nmz
	do i=1,nmr
	write(33,*)Er(i,j),Ez(i,j)
	enddo
enddo
close(33)
```

The units of the `1T1.T7` is [cm] and [V/cm]. Units are changed to [m] and [V/m] in the T source code.



| Parameter Name | Units | Type   | Default | Description                                                  |
| -------------- | ----- | ------ | ------- | ------------------------------------------------------------ |
| zedge          | m     | double | 0.0     | global position.                                             |
| L              | m     | double | 0.0     | length of the whole TWS structure linac with couplers excluded. L should be $n\times Lcav$ |
| freq           | Hz    | double | 2856e6  | RF frequency.                                                |
| phase          | deg   | double | 0.0     | initial phase.                                               |
| fileid         |       | int    | None    | `fileid=1`,  field file is `1T1.T7`.                         |

For static field, set `freq=0.0` and `phase=0.0`.






### WATCH

-2 element.

| Parameter Name | Units | Type   | Default | Description                                         |
| -------------- | ----- | ------ | ------- | --------------------------------------------------- |
| zedge          | m     | double | 0.0     | global position                                     |
| filename_id    |       | int    | 80      | fort.80, should < 100. But avoid using 40,50,60,70. |
| sample_freq    |       | int    | 1       | sample out freq.                                    |

