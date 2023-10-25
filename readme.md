This documents describes the usage of the lattice parser, which transform lte.impt => ImpactT.in.

The parser is developed by: Biaobin Li
Email: biaobin@ustc.edu.cn


# 最近加的功能

- 112 `EMfldCyl` element, one new field formats are added: static electric field (Poisson). The old field format is added back as `datafmt=imptold`.
- distribution 6, cylinder uniform, one could use `nbunch=10` to divide the bunch to 10 slices.



# 前期准备

## RFcoefext and RFcoeflcls

RFcoefext 用于处理非对称场分布，如行波加速管出、入口的coupler段，代码先将场分布进行镜像对称 (偶延拓)，然后对扩展后的对称分布的场进行傅里叶级数分解，取的系数越多，对原场分布拟合的就越准确：

```bash
$RFcoefext
 Please enter rfdata filename:
rfdata.in
 zhalf:    5.24636984
 How many Fourier coeficients you want?
39
 The RF data number is:         2002   10.4927397       5.24636984       5.24374796E-03
 input # of data points:
500
```

RFcoeflcls 则用于处理对称场分布，如$2\pi/3$ 模行波场 3个 cell 对应的场分布。



## 行波场的描述方法

Impact-T 对于行波场的处理，需要用户自己定义两个重叠的驻波场，从而形成行波场分布：
$$
E_{tr}(z,t) = \frac{1}{\sin(\beta_0 d)}[E_s(z)\cos(wt+\theta+\beta_0d-\pi/2)+E_s(z+d)\cos(wt+\theta+\pi/2)]
$$
对于 $2\pi/3$ 模行波场而言，$\beta_0d=2\pi/3$，所以处理步骤为：

- 获得左移一个cell 距离$d$ 的驻波场分布 $E_s(z+d)$，此步骤对应于 `RFcoeflcls` 需定义`shift length`：

	```bash
	$RFcoeflcls
	 Please give rfdata name:
	rfdata.in
	 zdata1:    5.2463999999999997E-002  0.15739000000000000
	 The RF data number is:          201  0.10492600000000001       0.10492699999999999        5.2463000000000000E-004
	 input shift length:
	```

- 在ImpactT.in 输入文件中，幅度值需乘以系数 $1/\sin(2\pi/3)$

- 两行叠加的驻波场分布，相位分别增加 30deg 和 90 deg



上面的幅值，如参照 Parmela 的输入文件设置，还需要考虑 $E_0T$ 带来的系数。



##  扫相

Impact-T 采用 $\cos()$ 函数习惯。扫相程序为 PhaseOpt.py 程序。



## Parmela 和 Impact-T RF场转化关系

Parmela uses $E_0$ or $E_0T$ to set the travelling or standing RF wave . However, Impact-T uses the $E_{max}$ to define the field. 
$$
E_0 = \frac{1}{L}\int |E_z(z)| dz, \\

T=\frac{1}{E_0 L}\int E_z(z)[\cos(kz)+i\sin(kz)]dz
$$
$E_0$ 指的是平均场强，对于常用的 SLAC $2\pi/3$  模，其 3cell 的场分布如下图所示。$E_0$ 值可积分得到，$E_{max}=1$ 作归一化时，$E_0$值约为0.78。

$T$ 实际上是 $T(k)$， 即基模($k=2\pi/3d$) 所占比例。由下图可见，基模 harm1 对应的幅值约为 0.744，故可得 $T=0.74/0.78\approx 0.95$.

对于Parmela 和 Impact-T 的场分布而言，即会有 $E_0$ 和 $E_0T$ 的转化关系。对于驻波而言：
$$
E_{impt}=E_{parm}/E_0,
$$
对于行波而言：
$$
E_{impt}=E_{parm}/(E_0 T).
$$


![image-20230506151845731](pics/image-20230506151845731.png)





# Control parameters

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



## Nbunch

For DC-GUN, since the bunch is really long, and large energy spread is induced. `Nbunch` should be set properly to accurately perform the Lorentz transformation from lab frame to beam frame.  

The `total_charge` in each `ImpactTj.in` refers to the individual slice charge, not the total charge of the whole bunch. One can use `fort.60` and `fort.70`  3rd col (current [A]) to check whether you have set the proper values.



## perdlen

Code was changed to use perdlen to cutoff the long tail particles. Only one RF bucket is considered if perdlen is set properly.

```fortran
+        !biaobin comments: Perdlen could be used to
+        ! cutoff particles outsize (-0.5zleng,0.5zleng)
+        ptrange(5) = zcent-0.5d0*zleng/Scxlt
+        ptrange(6) = zcent+0.5d0*zleng/Scxlt
+        !ptrange(5) = 0.0
+        !ptrange(6) = zleng/Scxlt
```





# Beam parameters

All beam section parameters in `lte.impt` are listed:

| Parameter Name    | Units | Type   | Default | Description                                                  |
| ----------------- | ----- | ------ | ------- | ------------------------------------------------------------ |
| mass              | eV    | double | 0.511e6 | mass of the particle.                                        |
| charge            |       | double | -1.0    | -1 for electron.                                             |
| distribution_type |       | int    | 2       | 6D gaussian distribution. See more options in Ji’s manual.   |
| Np                |       | int    | 1e3     | particle number.                                             |
| total_charge      | C     | double | 1e-9    | charge of the beam.                                          |
|                   |       |        |         |                                                              |
| emit_x            | m rad | double | 0.0     | emitance.                                                    |
| emit_nx           | m rad | double | 0.0     | normalized emittance.                                        |
| beta_x            | m     | double | 1.0     | twiss para.                                                  |
| alpha_x           |       | double | 0.0     | twiss para.                                                  |
| sigx              | m     | double | 0.0     | rms bunch size.                                              |
| sigpx             |       | double | 0.0     | rms value of $\gamma\beta_x/\gamma_0\beta_0$                 |
| dx                | m     | double | 0.0     | offset for x                                                 |
| emit_y            | m rad | double | 0.0     | emittance.                                                   |
| emit_ny           | m rad | double | 0.0     | normalized emittance.                                        |
| beta_y            | m     | double | 1.0     | twiss para.                                                  |
| alpha_y           |       | double | 0.0     | twiss para.                                                  |
| sigy              | m     | double | 0.0     | rms bunch size.                                              |
| sigpy             |       | double | 0.0     | rms value of $\gamma\beta_y/\gamma_0\beta_0$                 |
| dy                | m     | double | 0.0     | offset for y                                                 |
| emit_z            | m rad | double | 0.0     | emittance.                                                   |
| emit_nz           | m rad | double | 0.0     | normalized emittance.                                        |
| beta_z            | m     | double | 1.0     | twiss para.                                                  |
| alpha_z           |       | double | 0.0     | twiss para.                                                  |
| sigz              | m     | double | 0.0     | rms bunch length. $z=ct$, the actual value in `ImpactT.in` is $\beta_0 ct$, this is done in the python level. |
| sigpz             | eV    | double | 0.0     | rms value of $\gamma\beta_z/\gamma_0\beta_0$                 |
| dz                | m     | double | 0.0     | offset for z, $dz=cdt$. Value is transformed in the python level to $dz=\beta_0cdt$. |

Users could either use twiss parameters to define initial beam distribution, or use rms values. For $\sigma_{ij}\neq0$ cases, please use twiss-para. 

In the definition of python level: $z=ct$.

For ijk, like `distribution_type=112`, the `sigx,sigy` actually is beam radius `r`, and `sigz` is `Lbunch` full length, which in transverse direction is circle uniform, in longitudinal direction is flat-top, and z is in $[-Lbunch, 0]$ range. `zscale` is automatically given `1e-9` value in the python code as the manual said.



## Distribution

### 112, cylinder uniform

112 could be used to generate cylinder uniform distribution.

$\sigma_{i,j,k}$ are actually refers to:
$$
 sigx \rightarrow r\\
 sigz \rightarrow L_{bunch}
$$



### 6, cylinder uniform

I added it recently. The following transformation is done in the Fortran source code:
$$
\sigma_x=\sigma_y \\
r=2\times \sigma_x \\
L_{bunch}=2\sqrt{3}\sigma_z
$$
So in `lte.impt` input file,
$$
sigx \rightarrow beam~~radius~(r)\\
sigz\rightarrow bunch~~length~(L_{bunch})
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
| scale |  | double | 1.0 | When read in field in [Gauss], scale=1e-4. |

File format:

```matlab
%[cm, Gauss ]
r1 r2 nr
z1 z2 nz

k=1;
for j=1:nz+1
    for i=1:nr+1
        br(j,i)=B(k,1);
        bz(j,i)=B(k,2);
        k=k+1;
    end
end

```

For `fileid=3`, the B field file `1T3.T7`.

`1T3.T7out` would be generated, five columns data:

```bash
s[m], Br(r=0) [gauss], Br(r=0+dr) [gauss], Bz(r=0) [gauss], Bz(r=dr) [gauss]
```



Notes:

- The manual said V2 is not used, actually it is used as the `scale` value:

	```fortran
	!see SOl.f90/getfldt_sol()
	extfld(4) = scale*br*pos(1)/rr
	extfld(5) = scale*br*pos(2)/rr
	extfld(6) = scale*bz
	```

- If the field range in r-direction is not large enough, in case the particle is outside the B-field, program will stop:

	```fortran
	          !ir=r/hr+1, hr is dr 
	          if(ir.gt.fldata%NrIntvRft) then
	             print*,"ir: ",ir,rr,pos(1),pos(2),fldata%NzIntvRft,&
	             fldata%NrIntvRft,fldata%ZmaxRft,fldata%ZminRft,fldata%RmaxRft,&
	             fldata%RminRft
	!             print*,"ir: ",ir,rr,pos(1),pos(2)
	             stop
	          endif
	```

	

- Unit in IMPACT-T fortran code, [T] is used, Ji use `scale` to transform `Gauss` to `T`, when read-in file is [cm, gauss].



### SOLRF

105 element.

| Parameter Name | Units | Type   | Default | Description                   |
| -------------- | ----- | ------ | ------- | ----------------------------- |
| zedge          | m     | double | 0.0     | global position.              |
| L | m | double | 0.0 | the longitudinal length of the element (Blength). |
| Emax       | V/m | double | 1.0     | the absolute maximum values of on-axis Ez field. |
| freq           | Hz    | double | 2856e6  | frequency                     |
| phase          | deg   | double | 0.0     | RF design phase               |
| radius | m | double | 1.0 | pipe radius. |
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

~~The `z1,z2` and `length` in rfdatax will be updated according to `Lcell` and `Ncell`~~.



File format for `rfdatax`, unit is [m]:

```
39        /# of fourier coef. of Ez on axis, ncoefreal=20
z1        / [m]
z2        / z2-z1 could be N*Lfourier [m]
Lfourier  / one-period length [m]
a0
a1
b1
a2
b2
...
a19
b19
1.0       / in case B field is 0
0.0
0.0
0.0
0.0
```

In case only B field, like solenoid file:

```
1.0       / E field is 0
0.0
0.0
0.0
0.0
39        /# of fourier coef. of Bz on axis, ncoefreal=20
z0 
z1        / z1-z0 could be N*Lfourier
Lfourier  / one-period length
a0
a1
b1
a2
b2
...
a19
b19
```



### TWS

Traveling wave structure. Combined `SOLRF` elements.

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
| radius         | m     | double | 1.0      | pipe radius.                                                 |
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

112 element.

| Parameter Name | Units | Type   | Default | Description                                                  |
| -------------- | ----- | ------ | ------- | ------------------------------------------------------------ |
| zedge          | m     | double | 0.0     | global position.                                             |
| L              | m     | double | 0.0     | length of the whole TWS structure linac with couplers excluded. L should be $n\times Lcav$ |
| freq           | Hz    | double | 2856e6  | RF frequency.                                                |
| phase          | deg   | double | 0.0     | initial phase.                                               |
| fileid         |       | int    | None    | `fileid=1`,  field file is `1T1.T7`.                         |
| scale          |       | double | 1.0     | scale of the field strength.                                 |
| datafmt        |       | str    | “impt”  | the value could be `impt, imptold, poisson`. In `ImpactT.in`, corresponding values are `1,2,3`. The last one is Parmela’s static E-field format. |

The data format for `datafmt=impt,imptold,poisson` are listed as following:

```fortran
!datafmt=impt,  units are [cm,MHz,MV/m,A/m]
!--------------------------------------------
open(33,file=’1T1.T7’)
write(33,*)zmin,zmax,nmz-1
write(33,*)freq
write(33,*)rmin,rmax,nmr-1
do i=1,nmr
	do j=1,nmz
		write(33,*)Ez(j,i),Er(j,i),E(j,i),H(j,i)
	enddo
enddo
close(33)

!old version of IMPACT-T
!datafmt=imptold,  units are [cm,MHz,MV/m,A/m]
!--------------------------------------------
open(33,file=’1T1.T7’)
write(33,*)zmin,zmax,nmz-1
write(33,*)freq
write(33,*)rmin,rmax,nmr-1
do i=1,nmr
	do j=1,nmz
		write(33,*)Er(j,i),Ez(j,i),Etheta(j,i)
		write(33,*)Htheta
	enddo
enddo
close(33)

!=============================================
!datafmt=poisson, units are [cm,V/cm]
!-----------------------------------
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

Pay attention to that, for `impt,imptold`, it’s z-direction first, and then r-direction for data sampling. However, for `poisson`, which is Parmela static E-field or magnetic field format (`Poisson`), it’s r-direction first and then z-direction.

A new parameter `V12` is added for `112` element to determine the `datafmt`.



`1T1.T7out` would be generated for `datafmt=poisson`, three columns data:

```bash
s(m), Ez(r=0)[V/m], Er(r=rmax)[V/m] 
```



For electron DC gun:

```
dcgun: emfldcyl,zedge=0.0,L=6.5e-2,freq=0,phase=0,fileid=1, datafmt="poisson",scale=1.0
```

Parmela DC gun `ELEGUN1.T7`  field could be directly imported.



For 3D standing wave cavity, such as prebuncher:

```
!E0=0.19MV/m
!Emax=0.19/0.15=1.27MV/m

phipreb=93.46
Emaxfile=6.6479
Emaxpreb=0.19/0.150/Emaxfile
preb: emfldcyl, zedge=45e-2, L=20e-2, freq=476e6, phase=phipreb, fileid=3, scale=Emaxpreb, datafmt="impt"
```

Parmela prebuncher `CField` file could be directly imported.



### WATCH

-2 element.

-9 element.

| Parameter Name | Units | Type   | Default | Description                                         |
| -------------- | ----- | ------ | ------- | --------------------------------------------------- |
| zedge          | m     | double | 0.0     | global position                                     |
| filename_id    |       | int    | 80      | fort.80, should < 100. But avoid using 40,50,60,70. |
| sample_freq    |       | int    | 1       | sample out freq.                                    |
| slice_bin      |       | int    | 128     | slice bin for -9 element. See fort.180 file.        |



For slice information:

z(m), #, current(A), enx, eny, 

```Fortran
write(nfile,777)zz*scxlt,count(i),count(i)/(hz*scxlt)*sclcur*bet,epx(i)*scxlt,&
                epy(i)*scxlt,gam(i)*pmass,gam2uncor2(i)*pmass
```



### RCOL

-11 element, collimate particles with rectangular aperture.

| Parameter Name | Units | Type   | Default | Description                                        |
| -------------- | ----- | ------ | ------- | -------------------------------------------------- |
| zedge          | m     | double | 0.0     | global position                                    |
| xmax           | m     | double | 1.0     | xmax for x direction.                              |
| ymax           | m     | double | 1.0     | ymax for y direction.                              |
| xmin           | m     | double | None    | By default, xmin=None, then xmin=-xmax is applied. |
| ymin           | m     | double | None    | By default, ymin=None, then ymin=-ymax is applied. |



### MergeBin

-7 element

| Parameter Name | Units | Type   | Default | Description     |
| -------------- | ----- | ------ | ------- | --------------- |
| zedge          | m     | double | 0.0     | global position |



### SC2Dto3D

-5 element

| Parameter Name | Units | Type   | Default | Description     |
| -------------- | ----- | ------ | ------- | --------------- |
| zedge          | m     | double | 0.0     | global position |



### changedt

-4 element

| Parameter Name | Units | Type   | Default | Description     |
| -------------- | ----- | ------ | ------- | --------------- |
| zedge          | m     | double | 0.0     | global position |
| dt             | s     | double | 1e-12   | time step size  |



# Output

## fort.18

当前版本，对于孔径的处理，仍然只是全局孔径的设置：

<img src="pics/image-20231024102205713.png" alt="image-20231024102205713" style="zoom:67%;" />





## fort.40 and fort.50.







## slice info

`fort.60` and `fort.70`, the initial 3rd col. is $I/\beta_0$, that’s why for low energy beam, it’s very large. Now the source code is updated, the 3rd col is $I$ [A] now.

changes:

```fortran
-          write(nfile,777)zz*scxlt,count(i),count(i)/(hz*scxlt)*sclcur,epx(i)*scxlt,&
+          write(nfile,777)zz*scxlt,count(i),count(i)/(hz*scxlt)*sclcur*bet,epx(i)*scxlt,&
```







