# SANISand04
matlab code for SANISand04
* parameters define in this SANISand04
* note: this code can not simulation from p<=0kPa

```matlab
% undrained cyclic p=250,q=114.2,e_ini=0.737
% fine_sand_undrained p=100,1000,2000,3000kPa e_ini=0.735
% loose_sand_undrained p=100,1000,2000,3000kPa e_ini=0.833
% high_pressure_drained p=500kPa e_ini=0.810,0.886,0.960
% low_pressure_drained p=100kPa e_ini=0.831,0.917,0.996

G0       = 125.0;   
nu       = 0.05;
MC       = 1.25;
c_alpha  = 0.712;
lambda   = 0.019;
e0       = 0.938;
xi       = 0.7;
m        = 0.01;
h0       = 7.05;
c_h      = 0.968;
n_b      = 1.1;
A_0      = 0.704;
n_d      = 3.5;
z_max    = 5;
c_z      = 600.0;
p_a      = 101;
e_ini    = 0.735;
TolF     = 1.0e-7; % tolerence of yield function
tolR     = 1.0e-7; % tolerence of modified euler and RungeKutta45
small    = 1.0e-10;
nSub     = 10;
```

* reference：

  1. Dafalias Y F, Manzari M T. Simple plasticity sand model accounting for fabric change effects[J]. Journal of Engineering mechanics, 2004, 130(6): 622-634.
  2. Sloan S W, Abbo A J, Sheng D. Refined explicit integration of elastoplastic models with automatic error control[J]. Engineering Computations, 2001, 18(1/2): 121-194.

* special thanks to Alborz Ghofrani and Pedro Arduino for their open source code in OpenSees. 

  One can find their ManzariDafalias floder through this website. [OpenSees/OpenSees: OpenSees Source Code Repository (github.com)](https://github.com/OpenSees/OpenSees)

