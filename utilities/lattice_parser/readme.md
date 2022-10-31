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
| image_sc       |       | int    | 1       | Image charge flag.                                           |
| image_stop_pos | m     | double | 0.02    | position image charge forces are neglected.                  |
| meshx          |       | int    | 32      | number of mesh points in x direction.                        |
| meshy          |       | int    | 32      | number of mesh points in y direction.                        |
| meshz          |       | int    | 32      | number of mesh points in z direction.                        |
| Xrad           | m     | double | 0.015   | size of computational domain. Transverse size.               |
| Yrad           | m     | double | 0.015   | size of computational domain. Transverse size.               |
| PerdLen        | m     | double | 10.0    | PerdLen should be greater than the beam line lattice length. |
| Restart        |       | int    | 0       |                                                              |
| Nemission      |       | int    | 400     | the number of numerical emission steps.                      |
| Temission      | s     | double | 1e-9    | Laser pulse emission time.                                   |
| kinetic_energy | eV    | double | 0       | The real kinetic energy of the beam. NOT Bkenergy in line9.  |
| freq_rf_scale  | Hz    | double | 2856e6  | scale frequency $f_{scal}$,  $Scxl=c/(2\pi f_{scal}) $.      |
| ini_t          | s     | double | 0.0     | initial reference time.                                      |



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

Users could either use twiss parameters to define initial beam distribution, or use rms values. For $\sigma_{ij}\neq0$ cases, please use twiss-para. 



In the definition: $z=ct$.



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



### WATCH

-2 element.

| Parameter Name | Units | Type   | Default | Description      |
| -------------- | ----- | ------ | ------- | ---------------- |
| zedge          | m     | double | 0.0     | global position  |
| filename_id    |       | int    | 80      | fort.80          |
| sample_freq    |       | int    | 1       | sample out freq. |



### TWS

NOT added yet.

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
