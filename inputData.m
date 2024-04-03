%% ----------------------------SANISand04--------------------------------- %%
% Dafalias Y F, Manzari M T. Simple plasticity sand model accounting for 
% fabric change effects[J]. Journal of Engineering mechanics, 2004, 130(6): 622-634.
% https://doi.org/10.1061/(ASCE)0733-9399(2004)130:6(622)

%%
% 1. Input parameters for the model TO BE DEFINED
% data allocated in vector parms
%
% parms(1)     G0         Shear modulus constant
% parms(2)     nu         Poisson ratio
% parms(3)     Mc         CSL parameter, M_c
% parms(4)     c_alpha    CSL parameter, c=M_e/M_c
% parms(5)     lambda     CSL parameter
% parms(6)     e0         CSL parameter
% parms(7)     xi         CSL parameter
% parms(8)     m          Opening of yield surface cone
% parms(9)     h0         Plastic modulus constant
% parms(10)    c_h        Plastic modulus constant
% parms(11)    n_b        Plastic modulus constant
% parms(12)    A_0        Dilatancy parameter
% parms(13)    n_d        Dilatancy parameter
% parms(14)    z_max      Fabric-dilatancy parameter
% parms(15)    c_z        Fabric-dilatancy parameter
% parms(16)    p_a        Atmospheric pressure
% parms(17)    e_ini      Initial void ratio, (can modified)
% parms(18)    TolF       Tolerance for the yield function
% parms(19)    tolR       Tolerance for error in Runge-Kutta
% parms(20)    small      Small number
% parms(21)    nSub       Number of subincrements, use in intersectionFactor unloading 
% parms(22)    Pmin
% parms(23)    Presidual

%% parmeters given by paper
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
e_ini    = 0.737;
TolF     = 1.0e-7;
tolR     = 1.0e-7;
small    = 1.0e-10;
nSub     = 10;
Pmin     = 1.0e-4 * p_a;
Presidual = 1.0e-2 * p_a;
parms = [G0,nu,MC,c_alpha,lambda,e0,xi,m,h0,c_h,n_b,A_0,n_d,z_max,c_z,p_a,e_ini,TolF,tolR,small,nSub,Pmin,Presidual];

%% Another parameters, not used in this code, suit for SANISand-sf
% G0       = 125.0;   
% nu       = 0.05;
% MC       = 1.26;
% c_alpha  = 0.99;
% lambda   = 0.029;
% e0       = 0.78;
% xi       = 0.7;
% m        = 0.02;
% h0       = 5.0;
% c_h      = 0.968;
% n_b      = 1.15;
% A_0      = 0.35;
% n_d      = 1.25;
% z_max    = 18;
% c_z      = 500.0;
% p_a      = 101;
% e_ini    = 0.515;
% TolF     = 1.0e-7;
% tolR     = 1.0e-7;
% small    = 1.0e-10;
% nSub     = 10;
% Pmin     = 1.0e-4 * p_a;
% Presidual = 1.0e-2 * p_a;
% parms = [G0,nu,MC,c_alpha,lambda,e0,xi,m,h0,c_h,n_b,A_0,n_d,z_max,c_z,p_a,e_ini,TolF,tolR,small,nSub,Pmin,Presidual];

%%  recover variables for the model
m_P_atm = parms(16);
m_e_init = parms(17);
% mSig = [m_P_atm; m_P_atm; m_P_atm; 0.0; 0.0; 0.0]; % use this sometimes
mSig = [250; 250; 250; 0.0; 0.0; 0.0]; % initial stress

%%
% strain and stress terms and state variables
mEpsilon = zeros(6,1);    % Strain tensor, next(n+1) step
mEpsilon_n = zeros(6,1);  % Strain tensor at current(n) step,
mSigma = zeros(6,1);      % Stress tensor
mSigma_n = zeros(6,1);    % Stress tensor at current(n) step

mEpsilonE = zeros(6,1);   % Elastic strain
mEpsilonE_n = zeros(6,1); % Elastic strain, at current(n) step
mAlpha = zeros(6,1);      % Back stress ratio tensor, 
mAlpha_n = zeros(6,1);    % Back stress ratio tensor, at current(n) step
mAlpha_in = zeros(6,1);   % initial back stress ratio tensor at one load stage
mAlpha_in_n = zeros(6,1); % initial back stress ratio tensor at one load stage, at current(n) step
mDGamma = 0.0; 
mFabric = zeros(6,1);     % fabric tensor
mFabric_n = zeros(6,1);
mVoidRatio = m_e_init;

%% define
% NOTE: '_n' means at current step or now, 'm' infront the variables means
% material parameters

% 1. Stress tensor
% Stress tensor is allocated in matrix sig (1x6), voigt notation
% sig(1)  s11        11 -> 1
% sig(2)  s22        22 -> 2
% sig(3)  s33        33 -> 3
% sig(4)  s12        12 -> 4
% sig(5)  s23        23 -> 5 
% sig(6)  s13        13 -> 6

% 2. Strain tensor
% Strain tensor is allocated in matrix eps (1x6), voigt notation
% eps(1)  e11        11 -> 1
% eps(2)  e22        22 -> 2
% eps(3)  e33        33 -> 3
% eps(4)  2*e12      12 -> 4
% eps(5)  2*e23      23 -> 5
% eps(6)  2*e13      13 -> 6


% 3 State variables
% alpha(1)   alpha11     Back stress ratio tensor, alpha11
% alpha(2)   alpha22     Back stress ratio tensor, alpha22
% alpha(3)   alpha33     Back stress ratio tensor, alpha33
% alpha(4)   alpha12     Back stress ratio tensor, alpha12
% alpha(5)   alpha23     Back stress ratio tensor, alpha23
% alpha(6)   alpha13     Back stress ratio tensor, alpha13

% fabric(1)   fabric11    Fabric tensor, fabric11
% fabric(2)   fabric22    Fabric tensor, fabric22
% fabric(3)   fabric33    Fabric tensor, fabric33
% fabric(4)  fabric12    Fabric tensor, fabric12
% fabric(5)  fabric23    Fabric tensor, fabric23
% fabric(6)  fabric13    Fabric tensor, fabric13

% VoidRatio  void ratio   
