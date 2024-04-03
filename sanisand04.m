%% ----------------------------SANISand04--------------------------------- %%
% This is a Matlab implementation of the Sanisand 2004 model. The model is  
% developed by Dafalias and Manzari, implemented in the implicit integration 
% scheme, single Runge-Kutta method or modified Eluer(default) are not 
% accurate enough, MUST use stress correction.
% NOTE: 
% 1. The code is not fully validated, just show you some examples, use it at your own risk.
% 2. NOT work for p = 0 during the cyclic loading, need to be improved. See SANISand-sf for details.
% 3. The code is based on ManzariDafalias.cpp coded by Alborz Ghofrani, Pedro Arduino, just switch from C++ to Matlab.
%    One can find the original code from the website: https://github.com/OpenSees/OpenSees
% 4. The code is also inspired by the UMAT implementation and matlab code developed by C. Tamagnini, M. Martinelli, C. Miriano, 
%    modified by D. Masin. One can get these codes from the website: https://soilmodels.com/
% 5. The intergration scheme is based on paper: "Refined explicit integration of elastoplastic models with automatic error control" 
%    by Scott W. Sloan(2001), but not fully indentical.
% 6. Easy to understand, easy to use, easy to modify. One can easyly change to UMAT, VUMAT in ABAQUS.
% 7. The comments are mainly written by github copliot, not fully correct, but can help you understand the code.
% 8. The code is mainly written through VScode, so not obeys the MATLAB code style, but easy to read and understand with github copliot.
% 8. Written by Gl C, 2024.04.02, MATLAB R2022b, feels more convenient to use MATLAB than C++ and FORTRAN.
% -------------------------------------------------------------------------------------------------------------------------------------------------%

%%  Setup
clear all;
clc;
inputData;

%% calculate initial state
[mG, mK] = getElasticModuli(parms, mSig, m_e_init);  
mCe = getElasticStiffness(mG, mK);
mCep = mCe;
mCep_Consistent = mCe;
mUseElasticTan = false;
mElastFlag = 0; % 0:plastic, 1:elastic

%% Calculate stress/strain increment according to the state
% Monotonic drained
% Monotonic undrained
% cyclic undrained
%%


% %% Monotonic drained
% dStrain = zeros(6,1);
% Nstep = 250;
% epsa=-0.25;
% mSigma_n = mSig;
% servo_modulus = 1e7;
% sigma_t = mSig;
% dStrain = (sigma_t-mSigma_n)/servo_modulus; % acheive target stress by adjusting strain 
% dStrain(3) = -epsa/Nstep;
% sigma_tol = 0.001*getTrace(mSig); 
% % add initial state
% output1(1,:)=[mEpsilonE_n', mSigma_n', mAlpha_n', mFabric_n', mEpsilon_n', m_e_init, mElastFlag,0,0];
% IN = 1;
% while mEpsilon(3) <= -epsa/Nstep*250
%     fprintf('STEP: %d \n',IN);
%     mEpsilon = mEpsilon_n + dStrain; % +:compress
%     [mAlpha_in, mEpsilonE, mSigma, mAlpha, mFabric, mDGamma, mVoidRatio, mG, mK, mCe, mCep, mCep_Consistent,mElastFlag] =... 
%         intergration(parms,mSigma_n, mEpsilon_n, mEpsilonE_n, mAlpha_n, mFabric_n,mAlpha_in_n, mEpsilon, mElastFlag, mCe);
%     
%     mAlpha_in_n = mAlpha_in;
%     mSigma_n    = mSigma;
%     mEpsilon_n  = mEpsilon;
%     mEpsilonE_n = mEpsilonE;
%     mAlpha_n    = mAlpha;
%     mFabric_n   = mFabric;
%     sigma_t(3) = mSigma(3);
%     dSigma = sigma_t-mSigma;
%     dStrain = (sigma_t-mSigma)/servo_modulus;
%     if abs(dSigma(1)) < sigma_tol && abs(dSigma(2)) < sigma_tol
%         dStrain(3) = -epsa/Nstep;
%         output1(IN+1,:)=[mEpsilonE', mSigma', mAlpha', mFabric', mEpsilon',mVoidRatio,mElastFlag,0,0];
%         IN = IN + 1;
%     end
% end

%% Monotonic undrained
% dStrain = zeros(6,1);
% Nstep = 250;
% epsa=-0.25;
% dStrain(1)=0.5*epsa/Nstep;
% dStrain(2)=0.5*epsa/Nstep; 
% dStrain(3)=-epsa/Nstep; 
% dStrain(4)=0.0; 
% dStrain(5)=0.0;
% dStrain(6)=0.0;
% mSigma_n = mSig;
% 
% output1(1,:)=[mEpsilonE_n', mSigma_n', mAlpha_n', mFabric_n', mEpsilon_n', m_e_init, mElastFlag,0,0];
% for IN = 1:Nstep
%     fprintf('STEP: %d \n',IN);
%     mEpsilon = mEpsilon_n + dStrain; % +:compress
%     [mAlpha_in, mEpsilonE, mSigma, mAlpha, mFabric, mDGamma, mVoidRatio, mG, mK, mCe, mCep, mCep_Consistent,mElastFlag] =... 
%         intergration(parms,mSigma_n, mEpsilon_n, mEpsilonE_n, mAlpha_n, mFabric_n,mAlpha_in_n, mEpsilon, mElastFlag, mCe);
%     
%     mAlpha_in_n = mAlpha_in;
%     mSigma_n    = mSigma;
%     mEpsilon_n  = mEpsilon;
%     mEpsilonE_n = mEpsilonE;
%     mAlpha_n    = mAlpha;
%     mFabric_n   = mFabric;
%     output1(IN+1,:)=[mEpsilonE', mSigma', mAlpha', mFabric', mEpsilon',mVoidRatio,mElastFlag,0,0];
% end

%% cyclic undrained
dStrain = zeros(6,1);
Nstep = 2500;
epsa=-0.25;
dStrain(1)=0.5*epsa/Nstep;
dStrain(2)=0.5*epsa/Nstep; 
dStrain(3)=-epsa/Nstep; 
dStrain(4)=0.0; 
dStrain(5)=0.0;
dStrain(6)=0.0;
mSigma_n = mSig;

output1(1,:)=[mEpsilonE_n', mSigma_n', mAlpha_n', mFabric_n', mEpsilon_n', m_e_init, mElastFlag,0,1,mAlpha_in_n'];
cycle = 1;
IN = 1;
while cycle <= 40
    fprintf('STEP: %d, cycle: %d \n',IN,cycle);
    mEpsilon = mEpsilon_n + dStrain; % +:compress
    [mAlpha_in,mEpsilonE, mSigma, mAlpha, mFabric, mDGamma, mVoidRatio, mG, mK, mCe, mCep, mCep_Consistent,mElastFlag] =... 
        intergration(parms,mSigma_n, mEpsilon_n, mEpsilonE_n, mAlpha_n, mFabric_n, mAlpha_in_n, mEpsilon, mElastFlag, mCe);
    
    mAlpha_in_n = mAlpha_in;
    mSigma_n    = mSigma;
    mEpsilon_n  = mEpsilon;
    mEpsilonE_n = mEpsilonE;
    mAlpha_n    = mAlpha;
    mFabric_n   = mFabric;
    [p,q] = inv_s(mSigma);
    output1(IN+1,:)=[mEpsilonE', mSigma', mAlpha', mFabric', mEpsilon',mVoidRatio,mElastFlag,q,cycle,mAlpha_in'];
    if abs(q)>=114.2 % 114.2 65
        dStrain = -dStrain;
        cycle = cycle + 1;
    end

%     if mod(IN, 529)==0
%         pause;
%     end
%     if IN == 529
%         break;
%     end
    IN = IN +1;
end

%% plot 
draw_picture(parms,output1);


%% 
function alpha = getAlpha(parms,flag,n,gth,psi)
% Purpose: compute $\alpha$ tensor,
% Arguments:
%                           I/O   Type
%  parms                     I    R    : model constants given by the user
%  flag                      I    I    : flag to select the tensor to be computed
%  n                         I    6x1  : normal to the yield surface
%  gth                       I    R    : g(\theta,c) function
%  psi                       I    R    : state value parameter
%  alpha_c                   O    6x1  : flag==1, tensor $\alpha$ for the critical state cone
%  alpha_b                   O    6x1  : flag==2, tensor $\alpha$ for the bounding surface cone
%  alpha_d                   O    6x1  : flag==3, tensor $\alpha$ for the dilatancy cone

% Local variables
M   = parms(3);
m   = parms(8);
n_c = 0;
n_b = parms(11);
n_d = parms(13);
constant_n=[n_c,n_b,n_d];

% Select which alpha tensor to evaluate
if flag==1
      sgn=0; % sgn can be any value, since n_a=0
elseif flag==2
      sgn=-1;
elseif flag==3
      sgn=1;
else
    error('wrong flag, flag must be 1, 2 or 3');
end

alpha_th=gth*M*exp(sgn*constant_n(flag)*psi)-m;
alpha=sqrt(2/3) * alpha_th * n;
end


%%
function psi = getPSI(parms, void, p)
% Purpose: compute the state value parameter $\psi$ as a function of void ratio
% Arguments:
%                           I/O   Type
%  parms                     I    R    : model constants given by the user
%  void                      I    R    : current void ratio
%  p                         I    R    : current mean effective stress
%  psi                       O    R    : state value parameter

% Local variables
lambda = parms(5);
e0     = parms(6);
xi     = parms(7);
p_a    = parms(16);

% evaluate psi
ecrit = e0 - lambda*(p/p_a)^xi;
psi = void - ecrit;
end


%%
function cos3theta = getLodeAngle(n)
% Purpose: calculate cos(3*theta) from normal to the yield surface
% Arguments:
%                           I/O   Type
%  n                         I    6x1  : normal to the yield surface, contravariant form
%  cos3theta                 O    R    : cos(3*theta)

cos3theta = sqrt(6.0) * getTrace(singleDot(n,singleDot(n,n)));

if cos3theta > 1.0
    cos3theta=1.0;
elseif cos3theta < -1.0
    cos3theta=-1.0;
end
end


%%
function [gth,dgdth] = g_grad_g(parms,cos3theta)
% Purpose: calculate cos(3*theta) from deviatoric stress ratio tensor, not used
%   a) Argyris function (Argyris=.true.) (as in the original paper)
%      gth = two * c / ((one + c) - (one - c) * cos3theta)
%      dgdth = (one - c) * gth / (two * c)
%   b) Van Eekelen function (Argyris=.false.) (more appropriate for high friction angles)
%      gth = alpha / ((one + beta * cos3theta)^-0.25
%      dgdth = 1/g dgd(cos3teta) = -0.25 * beta * (one + beta * cos3theta)^-1

% Arguments:
%                           I/O   Type
%  cM                        I    R    : CSL parameter
%  cos3theta                 I    R    : cos(3*theta)
%  gth                       O    R    : g(\theta,c) function
%  dgdth                     O    R    : (1/g)dg/d\theta

% Local variables
Argyris = true;
cM = parms(4);
n_VE = -0.25d0;

% compute g function and its derivative
if Argyris    
      % Argyris function
      gth   = 2*cM/((1+cM)-(1-cM)*cos3theta);
      dgdth = (1-cM)*gth/(2*cM);
else  
      % Van Eekelen function
      alpha = (1+cM^(1/n_VE))^n_VE;
      beta  = (1-cM^(1/n_VE))/(1+cM^(1/n_VE));
      gth   = alpha*(1+beta*cos3theta)^(n_VE);
      dgdth = n_VE*alpha/(1+beta*cos3theta);
end
end


%%
function g = getG(cos3theta, c)
% Purpose: compute the g function
% Arguments:
%                           I/O   Type
%  cos3theta                 I    R    : cos(3*theta)
%  c                         I    R    : CSL parameter
%  g                         O    R    : g(\theta,c) function

% Local variables
g = 2 * c / ((1 + c) - (1 - c) * cos3theta);
end


%%
function yf = getF(parms,stress, alpha)
% Purpose: compute the yield function
% Arguments:
%                           I/O   Type
%  parms                     I    R    : model constants given by the user
%  stress                    I    6x1  : stress tensor, contravariant form
%  alpha                     I    6x1  : alpha tensor, contravariant form
%  yf                        O    R    : yield function
% 
% Local variables
m_m=parms(8); % opening of yield surface cone
m_Presidual = parms(23);
p=getTrace(stress)/3.0+m_Presidual;
s=getDevPart(stress);
temp = s - p*alpha; % s-p*alpha
% compute yield function
yf = getNormalContr(temp) - sqrt(2.0/3.0)*p*m_m;
end


%%
function dfdsig = getDfdsig(parms,stress, alpha)
% Purpose: compute the derivative of the yield function with respect to the
% stress tensor, not worked.
% Arguments:
%                           I/O   Type
%  parms                     I    R    : model constants given by the user
%  stress                    I    6x1  : stress tensor, contravariant form
%  alpha                     I    6x1  : alpha tensor, contravariant form
%  dfdsig                    O    6x1  : derivative of the yield function with respect to the stress tensor

% Local variables
p_limit_min = 1.0d-10;
n_limit_min = p_limit_min;
mI1=[1.0; 1.0; 1.0; 0.0; 0.0; 0.0];

% compute n, n:n = 1, but not a unit matrix when very low stress
stress_dev = getDevPart(stress);
tau=stress_dev-p*alpha;
norm_tau = (tau(:)'*tau(:))^0.5;
if norm_tau < n_limit_min
    norm_tau = n_limit_min;
end
n=tau/norm_tau;
% compute derivative of the yield function with respect to the stress tensor
p=getTrace(stress)/3.0;
if abs(p) < p_limit_min
    r = stress_dev/p_limit_min;
else
    r = stress_dev/p;
end
dfdsig = n - 1.0/3.0*(doubleDot2_2_Contr(n,r))*mI1;
end


%%
function [E_G, E_K]=getElasticModuli(parms, stress, void_ratio)
% Purpose: compute elastic moduli
% Arguments:
%                           I/O   Type
%  stress                    I    6x1  : stress tensor, contravariant form
%  parms                     I    R    : model constants given by the user
%  void_ratio                I    R    : void ratio

% Local variables
m_G0  = parms(1);
m_nu  = parms(2);
m_P_atm = parms(16);
m_Pmin = parms(21);
% compute mean effective stress
pn=getTrace(stress)/3.0;

if pn <= m_Pmin
    pn = m_Pmin;
end
% compute elastic moduli G and K
E_G = m_G0 * m_P_atm * (2.97 - void_ratio)^2.0 / (1.0 + void_ratio) * sqrt(pn / m_P_atm);
E_K = 2.0 / 3.0 * (1.0 + m_nu)/(1.0 - 2.0 * m_nu) * E_G;
end


%%
function De = getElasticStiffness(G, K)
% Purpose: compute elastic stiffness matrix
% Arguments:
%                           I/O   Type
%  G                         I    R    : shear modulus
%  K                         I    R    : bulk modulus
%  De                        O    6x6  : elastic stiffness matrix

% Local variables
F1  = K+4./3.*G;
F2  = K-2./3.*G;
De = [0. 0. 0. 0. 0. 0.; 0. 0. 0. 0. 0. 0.;0. 0. 0. 0. 0. 0.;...
      0. 0. 0. 0. 0. 0.;0. 0. 0. 0. 0. 0.;0. 0. 0. 0. 0. 0.]; 
% compute elastic stiffness matrix
for i = 1:1:3
    for j=1:1:3
        De(i,j) = F2;
    end
    De(i,i) = F1;
    De(i+3,i+3) = G;
end
end


%%
function D = getCompliance(G, K)
% Purpose: returns the compliance matrix in its covariant-covariant form
% Arguments:
%                           I/O   Type
%  G                         I    R    : shear modulus
%  K                         I    R    : bulk modulus
%  D                         O    6x6  : compliance matrix

% Local variables
D = zeros(6, 6);
a = 1 / (9*K) + 1 / (3*G);
b = 1 / (9*K) - 1 / (6*G);
c = 1 / G;

D(1,1) = a;
D(2,2) = a;
D(3,3) = a;
D(4,4) = c;
D(5,5) = c;
D(6,6) = c;
D(1,2) = b;
D(1,3) = b;
D(2,3) = b;
D(2,1) = b;
D(3,1) = b;
D(3,2) = b;
end


%%
function Dep = getElastoPlasticTangent(parms, nextStress, nextDGamma, G, K, B, C, D, h, n, b)     
% Purpose: compute elastoplastic tangent stiffness matrix
% Arguments:
%                           I/O   Type
%  parms                     I    R    : model constants given by the user
%  nextStress                I    6x1  : stress tensor
%  nextDGamma                I    R    : loading index, or plastic multiplier
%  G                         I    R    : shear modulus
%  K                         I    R    : bulk modulus
%  B                         I    R    : R parameter，B=1+\frac{3}{2}\frac{1-c}{c}g\cos3\theta 
%  C                         I    R    : R parameter，C=3\sqrt{\frac{3}{2}}\frac{1-c}{c}g
%  D                         I    R    : Dialtancy parameter, D=A_{d}(\mathbf{\alpha}_{\theta}^{d}-\mathbf{\alpha}):\mathbf{n}=A_{d}\mathbf{d}:\mathbf{n}
%  h                         I    R    : hardening parameter, h=\frac{b_0}{(\alpha-\alpha_{\mathrm{in}}):\mathbf{n}}
%  n                         I    6x1  : normal to the yield surface
%  b                         I    6x1  : distance from the boundary surface, \mathbf{b}=\mathbf{\alpha}_{\theta}^{b}-\mathbf{\alpha},
%  Dep                       O    6x6  : elastoplastic tangent stiffness matrix

% Local variables
one3 = 1.0/3.0;
small = parms(20);
m_Presidual = parms(23);
mI1=[1.0; 1.0; 1.0; 0.0; 0.0; 0.0];
p = one3 * getTrace(nextStress)+m_Presidual;
p = max(p, small + m_Presidual);
stress_dev = getDevPart(nextStress);
r = stress_dev/p;
Kp = 2.0/3.0*p*h*doubleDot2_2_Contr(b,n);
De = getElasticStiffness(K, G);  % De is the elastic stiffness matrix, aC
temp0 = B * n - C * (singleDot(n,n) - one3 * mI1) + one3 * D * mI1;
R=contravariant2Covariant(temp0);
temp1 = doubleDot4_2(De, R);

dfdsig = n-one3*doubleDot2_2_Contr(n,r)*mI1;
temp2 = contravariant2Covariant(dfdsig);
aTDe = temp2' * De;
temp3 = aTDe * R + Kp;
if abs(temp3) < small
    Dep = De;
else
    res1 = temp1 * aTDe;
    res2 = MacauleyIndex(nextDGamma);
    Dep = De - res2 * res1 / temp3;
end
end


%%
function [n, d, b, cos3Theta, h, psi, alphaBtheta, alphaDtheta,b0, A_d, D, B, C, R] = ...
    getStateDependent(parms, stress, alpha, fabric, e, alpha_in)
% Purpose: compute state dependent variables
% Arguments:
%                           I/O   Type
%  parms                     I    R    : model constants given by the user
%  stress                    I    6x1  : stress tensor, contravariant form
%  alpha                     I    6x1  : alpha tensor, contravariant form
%  alpha_in                  I    6x1  : initial alpha tensor, contravariant form
%  e                         I    R    : void ratio
%  fabric                    I    6x1  : fabric tensor, contravariant form
%  n                         O    6x1  : normal to the yield surface
%  b                         O    6x1  : distance from the boundary surface, \mathbf{b}=\mathbf{\alpha}_{\theta}^{b}-\mathbf{\alpha}
%  D                         O    R    : Dialtancy parameter, D=A_{d}(\mathbf{\alpha}_{\theta}^{d}-\mathbf{\alpha}):\mathbf{n}
%  B                         O    R    : R parameter，B=1+\frac{3}{2}\frac{1-c}{c}g\cos3\theta
%  C                         O    R    : R parameter，C=3\sqrt{\frac{3}{2}}\frac{1-c}{c}g
%  R                         O    6x1  : R tensor
 
% Local variables
one3 = 1.0/3.0;
root23 = sqrt(2.0/3.0);
mI1=[1.0; 1.0; 1.0; 0.0; 0.0; 0.0];
m_G0 = parms(1);
m_Mc = parms(3);
m_c  = parms(4);
m_m  = parms(8);
m_h0 = parms(9);
m_ch = parms(10);
m_nb = parms(11);
m_A0 = parms(12);
m_nd = parms(13);
m_P_atm = parms(16);
small = parms(20);
m_Presidual = parms(23);

% calculate normal to the yield surface, n 
p = one3 * getTrace(stress) + m_Presidual;
if p < small
    p = small;
end
n = getNormalToYield(parms, stress, alpha);

% calculate hardening parameter, h
temp0 = alpha - alpha_in;
AlphaAlphaInDotN = doubleDot2_2_Contr(temp0, n);
b0 = m_G0 * m_h0 * (1.0 - m_ch * e) / sqrt(p / m_P_atm);

if (abs(AlphaAlphaInDotN) < small)
    h = 1.0e10;
else
    h = b0 / AlphaAlphaInDotN;
end
h = min(1.0e7,b0 / AlphaAlphaInDotN);

% compute b, D
psi = getPSI(parms, e, p);
cos3Theta = getLodeAngle(n);
% alphaBtheta = getG(cos3Theta, m_c) * m_Mc * exp(-1.0 * m_nb * psi) - m_m;
% alphaDtheta = getG(cos3Theta, m_c) * m_Mc * exp(m_nd * psi) - m_m;
g= getG(cos3Theta, m_c);
alphaBtheta = getAlpha(parms,2,n,g,psi); % tensor, not scalar
alphaDtheta = getAlpha(parms,3,n,g,psi); % tensor, not scalar
b = alphaBtheta - alpha;
d = alphaDtheta - alpha;
A_d = m_A0 * (1 + Macauley(doubleDot2_2_Contr(fabric, n)));
D = A_d * doubleDot2_2_Contr(d, n);
if (p < 0.05 * m_P_atm)
    D_factor = 1.0 / (1.0 + (exp(7.6349 - 7.2713 * p)));
else
    D_factor = 1.0;
end    
D = D * D_factor;

% R = B * n - C * (singleDot(n,n) - one3 * mI1) + one3 * D * mI1;
B = 1.0 + 1.5 * (1 - m_c)/ m_c * getG(cos3Theta, m_c) * cos3Theta;
C = 3.0 * sqrt(1.5) * (1 - m_c)/ m_c * getG(cos3Theta, m_c);
temp0 = B * n - C * (singleDot(n,n) - one3 * mI1) + one3 * D * mI1;

% R=contravariant2Covariant(temp0);
R = temp0;
end 


%%
function n = getNormalToYield(parms, stress, alpha)
% compute the normal to the yield surface, n
% Arguments:
%                           I/O   Type
%  stress                    I    6x1  : stress tensor, contravariant form
%  alpha                     I    6x1  : alpha tensor, contravariant form
%  n                         O    6x1  : normal to the yield surface

% Local variables
one3 = 1.0/3.0;
small = parms(20);
m_Presidual = parms(23);

p = one3 * getTrace(stress) + m_Presidual;
if p < small
    n=zeros(6,1);
else
    stress_dev = getDevPart(stress);
    r = stress_dev/p;
    n = r - alpha;
    normN = getNormalContr(n);
    n = n/normN;
end
end


%%
function a = intersectionFactor(parms, curStress, curStrain, nestStrain,curAlpha, a0, a1)
% Purpose: compute the intersection factor
% Arguments:
%                           I/O   Type
%  curStress                 I    6x1  : current stress tensor, contravariant form
%  curStrain                 I    6x1  : current strain tensor, contravariant form
%  curAlpha                  I    6x1  : current alpha tensor, contravariant form
%  nestStrain                I    6x1  : next strain tensor, contravariant form
%  a0                        I    R    : initial value of the intersection factor
%  a1                        I    R    : final value of the intersection factor
%  a                         O    R    : intersection factor, 

% Local variables
a = a0;
m_e_init = parms(17);
m_TolF = parms(18);
small = parms(20);
m_Presidual = parms(23);
strainInc = nestStrain - curStrain;
void_ratio0 = m_e_init - (1 + m_e_init) * getTrace(curStrain + a0 * strainInc);
[G, K] = getElasticModuli(parms,curStress,void_ratio0);
De0 = getElasticStiffness(G, K);
dsigma0 = a0 * De0 * strainInc;
tempstress = curStress + dsigma0;
f0 = getF(parms,tempstress, curAlpha);

void_ratio1 = m_e_init - (1 + m_e_init) * getTrace(curStrain + a1 * strainInc);
[G, K] = getElasticModuli(parms,curStress, void_ratio1);
De1 = getElasticStiffness(G, K);
dsigma1 = a1 * De1 * strainInc;
f1 = getF(parms,curStress + dsigma1, curAlpha);
% fprintf('Calculate the elastic stress ratio between (%.3e, %.3e) \n',a0, a1);
for i = 1:1:10
    a = a1 - f1 * (a1 - a0) / (f1 - f0);
    dsigma = a * De1 * strainInc;
    f = getF(parms,curStress + dsigma, curAlpha);
%     fprintf('step %d, yf = %.4e, a = %.4e \n',i,f,a);
    if abs(f) < m_TolF
        break;
    end
    if f * f0 < 0 % f and f0 have different signs, f1 is plastic
        a1 = a;
        f1 = f;
    else
        f1 = f1 * f0 / (f0 + f);
        a0 = a;
        f0 = f;
    end
    if i == 10
        a = 0.0;
        disp('Can not find alpha!');
        break;
    end
    
end
if a > 1 - small
    a = 1.0;
end
if a < small
    a = 0.0;
end
end


%%
function a = intersectionFactor_Unloading(parms,curStress, curStrain, nextStrain, curAlpha)
% Purpose: compute the intersection factor for unloading
% Arguments:
%                           I/O   Type
%  curStress                 I    6x1  : current stress tensor, contravariant form
%  curStrain                 I    6x1  : current strain tensor, contravariant form
%  curAlpha                  I    6x1  : current alpha tensor, contravariant form
%  nestStrain                I    6x1  : next strain tensor, contravariant form
%  a                         O    R    : intersection factor,
% Local variables
a0 = 0.0;
a1 = 1.0;
m_e_init = parms(17);
m_TolF = parms(18);
small = parms(20);
nSub = parms(21);
m_Presidual = parms(23);
strainInc = nextStrain - curStrain;
void_ratio = m_e_init - (1 + m_e_init) * getTrace(curStrain);
[G, K] = getElasticModuli(parms, curStress, void_ratio);
dsigma = getElasticStiffness(G, K) * strainInc;
for i = 1:1:nSub
    da = (a1 - a0)/2.0;
    a = a1 - da;
    f = getF(parms,curStress + a * dsigma, curAlpha);
    if f > m_TolF  %  means f is positive, update upper bound, a1
        a1 = a;
    elseif f < -m_TolF  % means state is elastic, find lower bound, a0, break
        a0 = a;
        break;
    else            % means f is zero, find ultimate a
        return;
    end
    if i == nSub
        a = 0.0;
        disp('Didn''t find alpha - Unloading!');
        return;
    end
end
a = intersectionFactor(parms, curStress, curStrain, nextStrain, curAlpha, a0, a1);
end


%%
function [mAlpha_in, mEpsilonE, mSigma, mAlpha, mFabric, mDGamma, mVoidRatio, mG, mK, mCe, mCep, mCep_Consistent,mElastFlag] =... 
  intergration(parms,mSigma_n, mEpsilon_n, mEpsilonE_n, mAlpha_n, mFabric_n, mAlpha_in_n, mEpsilon, mElastFlag, mCe)
% Purpose: compute the response
% Arguments:
%                           I/O   Type
%  parms                     I    R    : model constants given by the user
%  curStress                 I    6x1  : current stress tensor, contravariant form
%  curStrain                 I    6x1  : current strain tensor, contravariant form
%  curAlpha                  I    6x1  : current alpha tensor, contravariant form
%  curFabric                 I    6x1  : current fabric tensor, contravariant form
%  nextStress                I    6x1  : next stress tensor, contravariant form
%  nextStrain                I    6x1  : next strain tensor, contravariant form
%  nextAlpha                 I    6x1  : next alpha tensor, contravariant form
%  nextFabric                I    6x1  : next fabric tensor, contravariant form

% Local variables
m_TolF = parms(18);
small = parms(20);
m_Presidual = parms(23);
mI1=[1.0; 1.0; 1.0; 0.0; 0.0; 0.0];

% update mAlpha_in, if unloadi 
temp0 = mEpsilon - mEpsilon_n;
trialDirection = mCe * temp0;
temp1 = mAlpha_n - mAlpha_in_n;
if doubleDot2_2_Contr(temp1, trialDirection) < 0.0
    mAlpha_in = mAlpha_n;
else
    mAlpha_in = mAlpha_in_n;
end
f_trial = getF(parms,trialDirection+mSigma_n, mAlpha_n);
if f_trial <= m_TolF
    mElastFlag = 1;
else
    mElastFlag = 0;
end

%  Force elastic response
if mElastFlag == 1
    [mEpsilonE, mSigma, mAlpha, mVoidRatio, mG, mK, mCe, mCep, mCep_Consistent ] = elastic_integrator(parms, mSigma_n, mEpsilon_n, ...
        mEpsilonE_n, mEpsilon);  
    mFabric = mFabric_n;
    mDGamma = 0;
else
    [mEpsilonE, mSigma, mAlpha, mFabric, mDGamma, mVoidRatio, mG, mK, mCe, mCep, mCep_Consistent] = ...
    explicit_integrator(parms, mSigma_n, mEpsilon_n, mEpsilonE_n,mAlpha_n, mFabric_n, mAlpha_in, mEpsilon);
end
end


%%
function [nextElasticStrain, nextStress, nextAlpha, nextVoidRatio, G, K, aC, aCep, aCep_Consistent ] = ...
elastic_integrator(parms, curStress, curStrain, curElasticStrain, nextStrain)                             
% Purpose: compute the elastic response
% Arguments:
%                           I/O   Type
%  parms                     I    R    : model constants given by the user
%  curStress                 I    6x1  : current stress tensor, contravariant form
%  curStrain                 I    6x1  : current strain tensor, contravariant form
%  nextStress                I    6x1  : next stress tensor, contravariant form
%  nextStrain                I    6x1  : next strain tensor, contravariant form
%  nextAlpha                 I    6x1  : next alpha tensor, contravariant form

% Local variables
m_e_init = parms(17);
small = parms(20);
m_Presidual = parms(23);

% calculate elastic response
dStrain = nextStrain - curStrain; 
nextVoidRatio = m_e_init - (1 + m_e_init) * getTrace(nextStrain);
nextElasticStrain = curElasticStrain + dStrain;
[G, K] = getElasticModuli(parms, curStress, nextVoidRatio);
aC = getElasticStiffness(G, K);
aCep = aC;
aCep_Consistent = aC;
nextStress = curStress + aC * dStrain;

% update state variables
p = 1.0/3.0 * getTrace(nextStress) + m_Presidual;
if p > small
    nextAlpha = getDevPart(nextStress) / p;
end
end


%%
function [nextElasticStrain, nextStress, nextAlpha, nextFabric, nextDGamma, nextVoidRatio, G, K, aC, aCep, aCep_Consistent ]= ...
explicit_integrator(parms, curStress, curStrain, curElasticStrain,curAlpha, curFabric,alpha_in, nextStrain)
% Purpose: compute the explicit response
% Arguments:
%                           I/O   Type
%  parms                     I    R    : model constants given by the user
%  curStress                 I    6x1  : current stress tensor, contravariant form
%  curStrain                 I    6x1  : current strain tensor, contravariant form
%  curElasticStrain          I    6x1  : current elastic strain tensor, contravariant form
%  curAlpha                  I    6x1  : current alpha tensor, contravariant form
%  curFabric                 I    6x1  : current fabric tensor, contravariant form
%  alpha_in                  I    6x1  : initial alpha tensor, contravariant form
%  nextStrain                I    6x1  : next strain tensor, contravariant form
%  nextElasticStrain         I    6x1  : next elastic strain tensor, contravariant form
%  nextStress                I    6x1  : next stress tensor, contravariant form
%  nextAlpha                 I    6x1  : next alpha tensor, contravariant form
%  nextFabric                I    6x1  : next fabric tensor, contravariant form
%  nextVoidRatio             I    R    : next void ratio

% Local variables
one3 = 1.0/3.0;
m_e_init = parms(17);
mTolF = parms(18);
small = parms(20);
nSub = parms(21);
m_Pmin = parms(22);
m_Presidual = parms(23);

p_tr_pos = true;
mI1=[1.0; 1.0; 1.0; 0.0; 0.0; 0.0];
nextVoidRatio = m_e_init - (1 + m_e_init) * getTrace(nextStrain);
dStrain = nextStrain - curStrain;
nextElasticStrain = curElasticStrain + dStrain;
[G, K] = getElasticModuli(parms, curStress, nextVoidRatio);
aC = getElasticStiffness(G, K);
dSigma = doubleDot4_2(aC, dStrain);
nextStress = curStress + dSigma;
f                      = getF(parms,nextStress, curAlpha);
p                      = one3 * getTrace(nextStress) + m_Presidual;

if p < m_Presidual
    p_tr_pos = false;
end 
if p_tr_pos && f <= mTolF
%   This is a pure elastic loading/unloading
    nextAlpha         = curAlpha;
    nextFabric        = curFabric;
    nextDGamma        = 0;
    aCep              = aC; 
    aCep_Consistent   = aC;
else
    fn = getF(parms,curStress, curAlpha);
    pn = one3 * getTrace(curStress) + m_Presidual;
    if pn < m_Presidual
        disp('p_n < 0, This should have not happened!');
        nextStress = m_Presidual * mI1;%m_Pmin * mI1;
        nextAlpha=zeros(6,1);
        nextFabric=curFabric;
        nextDGamma = 0;
        aCep              = aC; 
        aCep_Consistent   = aC;
        return;
    end
%     fprintf('normal mean stress is %.4e, yield function is %.4e \n',pn, fn);
   
    if (fn > mTolF)
        % This is an illegal stress state! This shouldn't happen.
        disp('stress state outside the yield surface!');
        disp('SANISAND04 : Encountered an illegal stress state!');
        fprintf('                  fn = %f\n', fn);
%         [nextElasticStrain, nextStress, nextAlpha, nextFabric, nextDGamma, nextVoidRatio, G, K, aC, aCep, aCep_Consistent] = ...
%         RungeKutta45(parms, curStress, curStrain, curElasticStrain,curAlpha, curFabric,alpha_in, nextStrain);
        [nextElasticStrain, nextStress, nextAlpha, nextFabric, nextDGamma, nextVoidRatio, G, K, aC, aCep, aCep_Consistent] = ...
        ModifiedEuler(parms, curStress, curStrain, curElasticStrain,curAlpha, curFabric,alpha_in, nextStrain);

    elseif fn < -mTolF
        disp('This is a transition from elastic to plastic');
        elasticRatio = intersectionFactor(parms,curStress, curStrain, nextStrain, curAlpha, 0.0, 1.0);
        % dSigma         = doubleDot4_2(aC, elasticRatio*(nextStrain - curStrain));
        dElasStrain = dStrain * elasticRatio;
        dSigma = doubleDot4_2(aC, dElasStrain);
        tempStress = curStress + dSigma;
        tempStrain = curStrain + dElasStrain;
        tempElasticStrain = curElasticStrain + dElasStrain;
%         [nextElasticStrain, nextStress, nextAlpha, nextFabric, nextDGamma, nextVoidRatio, G, K, aC, aCep, aCep_Consistent] = ...
%         RungeKutta45(parms, tempStress, tempStrain, tempElasticStrain, curAlpha, curFabric,alpha_in, nextStrain);
        [nextElasticStrain, nextStress, nextAlpha, nextFabric, nextDGamma, nextVoidRatio, G, K, aC, aCep, aCep_Consistent] = ...
        ModifiedEuler(parms, tempStress, tempStrain, tempElasticStrain, curAlpha, curFabric,alpha_in, nextStrain);
    elseif abs(fn) < mTolF
        dSigmaNormal = getNormalContr(dSigma);
        if dSigmaNormal == 0
            dSigmaNormal=1.0;
        end
        n = getNormalToYield(parms, curStress, curAlpha);
        if doubleDot2_2_Contr(n, dSigma)/dSigmaNormal > (-sqrt(mTolF))
            disp('This is a pure plastic step'); 
%             [nextElasticStrain, nextStress, nextAlpha, nextFabric, nextDGamma, nextVoidRatio, G, K, aC, aCep, aCep_Consistent] = ...
%             RungeKutta45(parms, curStress, curStrain, curElasticStrain,curAlpha, curFabric,alpha_in, nextStrain);
            [nextElasticStrain, nextStress, nextAlpha, nextFabric, nextDGamma, nextVoidRatio, G, K, aC, aCep, aCep_Consistent] = ...
            ModifiedEuler(parms, curStress, curStrain, curElasticStrain,curAlpha, curFabric,alpha_in, nextStrain);
        else
            disp('This is an elastic unloading followed by plastic loading') ;
            elasticRatio = intersectionFactor_Unloading(parms,curStress, curStrain, nextStrain, curAlpha);
            % dSigma         = doubleDot4_2(aC, elasticRatio*(nextStrain - curStrain));
            dElasStrain = dStrain * elasticRatio;
            dSigma = doubleDot4_2(aC, dElasStrain);
            tempStress = curStress + dSigma;
            tempStrain = curStrain + dElasStrain;
            tempElasticStrain = curElasticStrain + dElasStrain;
%             [nextElasticStrain, nextStress, nextAlpha, nextFabric, nextDGamma, nextVoidRatio, G, K, aC, aCep, aCep_Consistent] = ...
%             RungeKutta45(parms, tempStress, tempStrain, tempElasticStrain, curAlpha, curFabric,alpha_in, nextStrain);
            [nextElasticStrain, nextStress, nextAlpha, nextFabric, nextDGamma, nextVoidRatio, G, K, aC, aCep, aCep_Consistent] = ...
            ModifiedEuler(parms, tempStress, tempStrain, tempElasticStrain, curAlpha, curFabric,alpha_in, nextStrain);
        end
    end
end  
end


%%
function [nextElasticStrain, nextStress, nextAlpha, nextFabric, nextDGamma, nextVoidRatio, G, K, aC, aCep, aCep_Consistent] = ...
  RungeKutta45(parms, curStress, curStrain, curElasticStrain, curAlpha, curFabric, alpha_in, nextStrain)    
% Purpose: RungeKutta45 scheme solve the elastoplastic problem
% Arguments:
%                           I/O   Type
%  parms                     I    R    : model constants given by the user
%  curStress                 I    6x1  : current stress tensor, contravariant form

% Local variables
T = 0.0;
dT = 1.0;
dT_min = 1.0e-3;
one3 = 1.0/3.0;
two3 = 2.0/3.0;
m_Mc = parms(3);
m_z_max = parms(14);
m_cz = parms(15);
m_e_init = parms(17);
TolE = parms(19);
small = parms(20);
m_Pmin = parms(22);
m_Presidual = parms(23);
mI1 = [1.0; 1.0; 1.0; 0.0; 0.0; 0.0];
%  4th order mixed variant identity tensor
mIImix = eye(6) * 1.0;
% Set aCep_Consistent to zero for substepping process
aCep = zeros(6,6);
aCep_Consistent = zeros(6,6);

curVoidRatio      = m_e_init - (1 + m_e_init) * getTrace(curStrain);
nextVoidRatio     = m_e_init - (1 + m_e_init) * getTrace(nextStrain);
nextElasticStrain = curElasticStrain + (nextStrain - curStrain);
[G, K] = getElasticModuli(parms, curStress, curVoidRatio);
aC = getElasticStiffness(G, K);
aD = getCompliance(G, K);

nextStress = curStress;
nextAlpha = curAlpha;
nextFabric = curFabric;
p = one3 * getTrace(nextStress);
if (p < m_Pmin)
    nextStress = getDevPart(nextStress) + m_Pmin * mI1;
    p = one3 * getTrace(nextStress);
end

% substepping process, record the step
iteration = 1;
while (T < 1.0)
    nextVoidRatio     = m_e_init - (1 + m_e_init) * getTrace(nextStrain + T * (nextStrain - curStrain));
    dVolStrain = dT * getTrace(nextStrain - curStrain);
    dDevStrain = dT * getDevPart(nextStrain - curStrain);
    % Calc Delta 1    
    thisSigma = nextStress;
    thisAlpha = nextAlpha;
    thisFabric = nextFabric;
    [n, d, b, cos3Theta, h, psi, alphaBtheta, alphaDtheta,b0, A_d, D, B, C, R] = getStateDependent(...
    parms, thisSigma, thisAlpha, thisFabric, nextVoidRatio, alpha_in);
    [G, K] = getElasticModuli(parms, thisSigma, nextVoidRatio);
    r = getDevPart(nextStress) / p;
    Kp = two3 * p * h * doubleDot2_2_Contr(b, n);
    temp4 = (Kp + 2.0*G*(B-C*getTrace(singleDot(n,singleDot(n,n)))) - K*D*doubleDot2_2_Contr(n,r));  
    if (abs(temp4) < small) 
        temp4 = small;
    end    
    nextDGamma  = (2.0*G*doubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*doubleDot2_2_Contr(n,r))/temp4;
    dSigma1     = 2.0*G* covariant2Contraviant(dDevStrain) + K*dVolStrain*mI1 - Macauley(nextDGamma)*...
                 (2.0*G*(B*n-C*(singleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
    dAlpha1     = Macauley(nextDGamma) * two3 * h * b;
    dFabric1    = -1.0 * Macauley(nextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + thisFabric);
    dPStrain1   = nextDGamma * contravariant2Covariant(R); % 
    aCep1 = getElastoPlasticTangent(parms, thisSigma, nextDGamma, G, K, B, C, D, h, n, b);

    % Calc Delta 2
    thisSigma = nextStress  + 0.5*dSigma1;
    thisAlpha = nextAlpha   + 0.5*dAlpha1;
    thisFabric = nextFabric + 0.5*dFabric1;
    [n, d, b, cos3Theta, h, psi, alphaBtheta, alphaDtheta,b0, A_d, D, B, C, R] = getStateDependent(...
    parms, thisSigma, thisAlpha, thisFabric, nextVoidRatio, alpha_in);
    [G, K] = getElasticModuli(parms, thisSigma, nextVoidRatio);
    r = getDevPart(nextStress) / p;
    Kp = two3 * p * h * doubleDot2_2_Contr(b, n);
    temp4 = (Kp + 2.0*G*(B-C*getTrace(singleDot(n,singleDot(n,n)))) - K*D*doubleDot2_2_Contr(n,r));
    if (abs(temp4) < small) 
        temp4 = small;
    end
    nextDGamma  = (2.0*G*doubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*doubleDot2_2_Contr(n,r))/temp4;
    dSigma2     = 2.0*G* covariant2Contraviant(dDevStrain) + K*dVolStrain*mI1 - Macauley(nextDGamma)*...
                (2.0*G*(B*n-C*(singleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
    dAlpha2     = Macauley(nextDGamma) * two3 * h * b;
    dFabric2    = -1.0 * Macauley(nextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + thisFabric);
    dPStrain2   = nextDGamma * contravariant2Covariant(R); % 
    aCep2 = getElastoPlasticTangent(parms, thisSigma, nextDGamma, G, K, B, C, D, h, n, b);

    % Calc Delta 3
    thisSigma = nextStress  + 0.25*(dSigma1  + dSigma2);
    thisAlpha = nextAlpha   + 0.25*(dAlpha1  + dAlpha2);
    thisFabric = nextFabric + 0.25*(dFabric1 + dFabric2);
    [n, d, b, cos3Theta, h, psi, alphaBtheta, alphaDtheta,b0, A_d, D, B, C, R] = getStateDependent(...
    parms, thisSigma, thisAlpha, thisFabric, nextVoidRatio, alpha_in);
    [G, K] = getElasticModuli(parms, thisSigma, nextVoidRatio);
    r = getDevPart(nextStress) / p;
    Kp = two3 * p * h * doubleDot2_2_Contr(b, n);
    temp4 = (Kp + 2.0*G*(B-C*getTrace(singleDot(n,singleDot(n,n)))) - K*D*doubleDot2_2_Contr(n,r));
    if (abs(temp4) < small) 
        temp4 = small;
    end
    nextDGamma  = (2.0*G*doubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*doubleDot2_2_Contr(n,r))/temp4;
    dSigma3     = 2.0*G* covariant2Contraviant(dDevStrain) + K*dVolStrain*mI1 - Macauley(nextDGamma)*...
                (2.0*G*(B*n-C*(singleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
    dAlpha3     = Macauley(nextDGamma) * two3 * h * b;
    dFabric3    = -1.0 * Macauley(nextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + thisFabric);
    dPStrain3   = nextDGamma * contravariant2Covariant(R); % 
    aCep3 = getElastoPlasticTangent(parms, thisSigma, nextDGamma, G, K, B, C, D, h, n, b);

    % Calc Delta 4
    thisSigma = nextStress  - dSigma2  + 2*dSigma3;
    thisAlpha = nextAlpha   - dAlpha2  + 2*dAlpha3;
    thisFabric = nextFabric - dFabric2 + 2*dFabric3;
    [n, d, b, cos3Theta, h, psi, alphaBtheta, alphaDtheta,b0, A_d, D, B, C, R] = getStateDependent(...
    parms, thisSigma, thisAlpha, thisFabric, nextVoidRatio, alpha_in);
    [G, K] = getElasticModuli(parms, thisSigma, nextVoidRatio);
    r = getDevPart(nextStress) / p;
    Kp = two3 * p * h * doubleDot2_2_Contr(b, n);
    temp4 = (Kp + 2.0*G*(B-C*getTrace(singleDot(n,singleDot(n,n)))) - K*D*doubleDot2_2_Contr(n,r));
    if (abs(temp4) < small) 
        temp4 = small;
    end
    nextDGamma  = (2.0*G*doubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*doubleDot2_2_Contr(n,r))/temp4;
    dSigma4     = 2.0*G* covariant2Contraviant(dDevStrain) + K*dVolStrain*mI1 - Macauley(nextDGamma)*...
                (2.0*G*(B*n-C*(singleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
    dAlpha4     = Macauley(nextDGamma) * two3 * h * b;
    dFabric4    = -1.0 * Macauley(nextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + thisFabric);
    dPStrain4   = nextDGamma * contravariant2Covariant(R); % 
    aCep4 = getElastoPlasticTangent(parms, thisSigma, nextDGamma, G, K, B, C, D, h, n, b);

    % Calc Delta 5
    thisSigma = nextStress  + (7*dSigma1  + 10*dSigma2  + dSigma4)/27;
    thisAlpha = nextAlpha   + (7*dAlpha1  + 10*dAlpha2  + dAlpha4)/27;
    thisFabric = nextFabric + (7*dFabric1 + 10*dFabric2 + dFabric4)/27;
    [n, d, b, cos3Theta, h, psi, alphaBtheta, alphaDtheta,b0, A_d, D, B, C, R] = getStateDependent(...
    parms, thisSigma, thisAlpha, thisFabric, nextVoidRatio, alpha_in);
    [G, K] = getElasticModuli(parms, thisSigma, nextVoidRatio);
    r = getDevPart(nextStress) / p;
    Kp = two3 * p * h * doubleDot2_2_Contr(b, n);
    temp4 = (Kp + 2.0*G*(B-C*getTrace(singleDot(n,singleDot(n,n)))) - K*D*doubleDot2_2_Contr(n,r));
    if (abs(temp4) < small) 
        temp4 = small;
    end
    nextDGamma  = (2.0*G*doubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*doubleDot2_2_Contr(n,r))/temp4;
    dSigma5     = 2.0*G* covariant2Contraviant(dDevStrain) + K*dVolStrain*mI1 - Macauley(nextDGamma)*...
                (2.0*G*(B*n-C*(singleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
    dAlpha5     = Macauley(nextDGamma) * two3 * h * b;
    dFabric5    = -1.0 * Macauley(nextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + thisFabric);
    dPStrain5   = nextDGamma * contravariant2Covariant(R); % 
    aCep5 = getElastoPlasticTangent(parms, thisSigma, nextDGamma, G, K, B, C, D, h, n, b);

    % Calc Delta 6
    thisSigma = nextStress  + (28*dSigma1  - 125*dSigma2  + 546*dSigma3  + 54*dSigma4  - 378*dSigma5)/625;
    thisAlpha = nextAlpha   + (28*dAlpha1  - 125*dAlpha2  + 546*dAlpha3  + 54*dAlpha4  - 378*dAlpha5)/625;
    thisFabric = nextFabric + (28*dFabric1 - 125*dFabric2 + 546*dFabric3 + 54*dFabric4 - 378*dFabric5)/625;
    [n, d, b, cos3Theta, h, psi, alphaBtheta, alphaDtheta,b0, A_d, D, B, C, R] = getStateDependent(...
    parms, thisSigma, thisAlpha, thisFabric, nextVoidRatio, alpha_in);
    [G, K] = getElasticModuli(parms, thisSigma, nextVoidRatio);
    r = getDevPart(nextStress) / p;
    Kp = two3 * p * h * doubleDot2_2_Contr(b, n);
    temp4 = (Kp + 2.0*G*(B-C*getTrace(singleDot(n,singleDot(n,n)))) - K*D*doubleDot2_2_Contr(n,r));
    if (abs(temp4) < small) 
        temp4 = small;
    end
    nextDGamma  = (2.0*G*doubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*doubleDot2_2_Contr(n,r))/temp4;
    dSigma6     = 2.0*G* covariant2Contraviant(dDevStrain) + K*dVolStrain*mI1 - Macauley(nextDGamma)*...
                (2.0*G*(B*n-C*(singleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
    dAlpha6     = Macauley(nextDGamma) * two3 * h * b;
    dFabric6    = -1.0 * Macauley(nextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + thisFabric);
    dPStrain6   = nextDGamma * contravariant2Covariant(R); % 
    aCep6 = getElastoPlasticTangent(parms, thisSigma, nextDGamma, G, K, B, C, D, h, n, b);

    % Update
    dSigma =   ( 14*dSigma1   + 35*dSigma4   + 162*dSigma5   + 125*dSigma6   ) / 336;
    dAlpha =   ( 14*dAlpha1   + 35*dAlpha4   + 162*dAlpha5   + 125*dAlpha6   ) / 336;
    dFabric =  ( 14*dFabric1  + 35*dFabric4  + 162*dFabric5  + 125*dFabric6  ) / 336;
    dPStrain = ( 14*dPStrain1 + 35*dPStrain4 + 162*dPStrain5 + 125*dPStrain6 ) / 336;

    nStress = nextStress + dSigma;
    nAlpha  = nextAlpha  + dAlpha;
    nFabric = nextFabric  + dFabric;
  
    % compute the error
    p = one3 * getTrace(nStress);
    if (p < 0)
        if (dT == dT_min)
            return;
        end
        dT = max(0.1 * dT, dT_min);
        continue;
    end
    stressNorm = getNormalContr(nextStress);
    alphaNorm = getNormalContr(nextAlpha);
    curStepError1 = getNormalContr(-42*dSigma1 - 224*dSigma3 - 21*dSigma4 + 162*dSigma5 + 125*dSigma6 )/336;
    if (stressNorm >= 0.5)
        curStepError1 = curStepError1 / (2 * stressNorm);
    end
    curStepError2 = getNormalContr(-42*dAlpha1 - 224*dAlpha3 - 21*dAlpha4 + 162*dAlpha5 + 125*dAlpha6 )/336;
    if (alphaNorm >= 0.5)
        curStepError2 = curStepError2 / (2 * alphaNorm);
    end
    curStepError = max(curStepError1, curStepError2);
   
    fprintf('iteration: %d, error: %.5e, T: %.5e, dT: %.5e \n',iteration, curStepError, T, dT);

    if (curStepError > TolE)
%         q = max(0.8 * sqrt(TolE / curStepError), 0.1);
        q = max(0.9 * (TolE / curStepError)^ 0.2, 0.1);
        fprintf('--- unsuccessful increment\n ');
        if (dT == dT_min)
            % mUseElasticTan = true;
            disp('Warning! Convergence not achieved until dT equal to DT_min.');
            nextElasticStrain = nextElasticStrain - dPStrain;
            nextStress = curStress;
            eta = sqrt(13.5) * getNormalContr(getDevPart(nextStress)) / getTrace(nextStress);
            if (eta > m_Mc)
                nextStress = one3 * getTrace(nextStress) * mI1 + m_Mc / eta * getDevPart(nextStress);
            end
            nextAlpha  = curAlpha + 3.0 * (getDevPart(nextStress)/getTrace(nextStress) - getDevPart(curStress)/getTrace(curStress));
            T = T + dT;
        end
        dT = max(q * dT, dT_min);
    else
        fprintf('+++ successful increment\n ');
        nextElasticStrain = nextElasticStrain - dPStrain;
        nextStress = nStress;
        nextAlpha  = nAlpha;
        nextFabric = nFabric;
        aCep_thisStep =   ( 14*aCep1 + 35*aCep4 + 162*aCep5 + 125*aCep6 ) /336;
        aCep_Consistent = aCep_thisStep * (aD * aCep_Consistent + T * mIImix);
        T = T + dT;
        [nextStress,nextAlpha,nextElasticStrain,aCep,aCep_Consistent] = stressCorrection(parms, curStress, curAlpha, curElasticStrain,...
        nextStress, nextAlpha, nextFabric,nextElasticStrain, nextVoidRatio, alpha_in, G,K,aCep, aCep_Consistent);
        if (curStepError == 0)
            dT = 1 - T;
        else
            q = max(0.9 * (TolE / curStepError)^ 0.2, 0.1);
            dT = max(q * dT, dT_min);
            dT = min(dT, 1 - T);
        end
    end
    iteration = iteration + 1;
end
fn = getF(parms,nextStress,nextAlpha) % 0.0999 mark
% aCep = aCep_thisStep; 
end


%%
function [nextElasticStrain, nextStress, nextAlpha, nextFabric, nextDGamma, nextVoidRatio, G, K, aC, aCep, aCep_Consistent] = ...
    ModifiedEuler(parms, curStress, curStrain, curElasticStrain, curAlpha, curFabric, alpha_in, nextStrain)
% Purpose: ModifiedEuler scheme solve the elastoplastic problem
% Arguments:
%                           I/O   Type
%  parms                     I    R    : model constants given by the user
%  curStress                 I    6x1  : current stress tensor, contravariant form
%  curStrain                 I    6x1  : current strain tensor, contravariant form
%  curElasticStrain          I    6x1  : current elastic strain tensor, contravariant form
%  curAlpha                  I    6x1  : current back stress tensor, contravariant form
%  curFabric                 I    6x1  : current fabric tensor, contravariant form
%  alpha_in                  I    R    : initial fabric anisotropy
%  nextStrain                I    6x1  : next strain tensor, contravariant form
%  nextElasticStrain         O    6x1  : next elastic strain tensor, contravariant form
%  nextStress                O    6x1  : next stress tensor, contravariant form
%  nextAlpha                 O    6x1  : next back stress tensor, contravariant form
%  nextFabric                O    6x1  : next fabric tensor, contravariant form
% Local variables
T = 0.0;
dT = 1.0;
dT_min = 1.0e-6; % mark
one3 = 1.0/3.0;
two3 = 2.0/3.0;
m_Mc = parms(3);
m_z_max = parms(14);
m_cz = parms(15);
m_e_init = parms(17);
TolE = 1e-4; % TolE = parms(19); mark
q = 1.0;
small = parms(20);
m_Pmin = parms(22);
m_Presidual = parms(23);
mI1 = [1.0; 1.0; 1.0; 0.0; 0.0; 0.0];
%  4th order mixed variant identity tensor
mIImix = eye(6) * 1.0;
% Set aCep_Consistent to zero for substepping process  

aCep = zeros(6,6);
aCep1 = zeros(6,6);
aCep2 = zeros(6,6);
aCep_Consistent = zeros(6,6);

curVoidRatio      = m_e_init - (1 + m_e_init) * getTrace(curStrain);
nextVoidRatio     = m_e_init - (1 + m_e_init) * getTrace(nextStrain);
nextElasticStrain = curElasticStrain + (nextStrain - curStrain);
[G, K] = getElasticModuli(parms, curStress, curVoidRatio);
aC = getElasticStiffness(G, K);
aD = getCompliance(G, K);

nextStress = curStress;
nextAlpha = curAlpha;
nextFabric = curFabric;
p = one3 * getTrace(nextStress);
if (p < m_Pmin + m_Presidual)
    nextStress = getDevPart(nextStress) + m_Pmin * mI1;
    p = m_Pmin;
end
iteration = 1;
while (T < 1.0)
    nextVoidRatio     = m_e_init - (1 + m_e_init) * getTrace(nextStrain + T * (nextStrain - curStrain));
    dVolStrain = dT * getTrace(nextStrain - curStrain);
    dDevStrain = dT * getDevPart(nextStrain - curStrain);
    % Calc Delta 1    
    thisSigma = nextStress;
    thisAlpha = nextAlpha;
    thisFabric = nextFabric;
    [n, d, b, cos3Theta, h, psi, alphaBtheta, alphaDtheta,b0, A_d, D, B, C, R] = getStateDependent(...
    parms, thisSigma, thisAlpha, thisFabric, nextVoidRatio, alpha_in);
    [G, K] = getElasticModuli(parms, thisSigma, nextVoidRatio);
    r = getDevPart(nextStress) / p;
    Kp = two3 * p * h * doubleDot2_2_Contr(b, n);
    temp4 = (Kp + 2.0*G*(B-C*getTrace(singleDot(n,singleDot(n,n)))) - K*D*doubleDot2_2_Contr(n,r));  
    if (abs(temp4) < small) 
        dSigma1 = zeros(6,1);
        dAlpha1 = zeros(6,1);
        dFabric1 = zeros(6,1);
        dPStrain1 = dDevStrain + dVolStrain*mI1;
    else    
        nextDGamma  = (2.0*G*doubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*doubleDot2_2_Contr(n,r))/temp4;
        if nextDGamma < 0.0
            nextDGamma = 0.0;
            dSigma1 = 2.0*G* covariant2Contraviant(dDevStrain) + K*dVolStrain*mI1;
            dAlpha1 = 3.0*(getDevPart(nextStress + dSigma1) / getTrace(nextStress + dSigma1) - getDevPart(nextStress) / getTrace(nextStress));
            dAlpha1 = zeros(6,1);
            dFabric1 = zeros(6,1);
            dPStrain1 = zeros(6,1);
        else
            dSigma1   = 2.0*G* covariant2Contraviant(dDevStrain) + K*dVolStrain*mI1 - Macauley(nextDGamma)*...
                            (2.0*G*(B*n-C*(singleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
            dAlpha1     = Macauley(nextDGamma) * two3 * h * b;
            dFabric1    = -1.0 * Macauley(nextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + thisFabric);
            dPStrain1   = nextDGamma * contravariant2Covariant(R);
        end 
        aCep1 = getElastoPlasticTangent(parms, thisSigma, nextDGamma, G, K, B, C, D, h, n, b);% marker
    end

    % Calc Delta 2
    thisSigma = nextStress  + dSigma1;
    thisAlpha = nextAlpha   + dAlpha1;
    thisFabric = nextFabric + dFabric1;
    p = one3 * getTrace(nextStress + dSigma1) + m_Presidual;
    if p < m_Presidual
        if dT == dT_min
            return;
        end
        dT = max(0.1 * dT, dT_min);
        continue;
    end
    [n, d, b, cos3Theta, h, psi, alphaBtheta, alphaDtheta,b0, A_d, D, B, C, R] = getStateDependent(...
    parms, thisSigma, thisAlpha, thisFabric, nextVoidRatio, alpha_in);
    [G, K] = getElasticModuli(parms, thisSigma, nextVoidRatio);
    r = getDevPart(nextStress) / p;
    Kp = two3 * p * h * doubleDot2_2_Contr(b, n);
    temp4 = (Kp + 2.0*G*(B-C*getTrace(singleDot(n,singleDot(n,n)))) - K*D*doubleDot2_2_Contr(n,r));
    if (abs(temp4) < small) 
        % neutral loading
        dSigma2 = zeros(6,1);
        dAlpha2 = zeros(6,1);
        dFabric2 = zeros(6,1);
        dPStrain2 = dDevStrain + dVolStrain*mI1;
    else
        nextDGamma  = (2.0*G*doubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*doubleDot2_2_Contr(n,r))/temp4;
        if nextDGamma < 0.0
            nextDGamma = 0.0;
            dSigma2 = 2.0*G* covariant2Contraviant(dDevStrain) + K*dVolStrain*mI1;
            dAlpha2 = 3.0*(getDevPart(nextStress + dSigma2) / getTrace(nextStress + dSigma2) - getDevPart(nextStress) / getTrace(nextStress));
            dFabric2 = zeros(6,1);
            dPStrain2 = zeros(6,1);
        else
            dSigma2   = 2.0*G* covariant2Contraviant(dDevStrain) + K*dVolStrain*mI1 - Macauley(nextDGamma)*...
                            (2.0*G*(B*n-C*(singleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
            dAlpha2     = Macauley(nextDGamma) * two3 * h * b;
            dFabric2    = -1.0 * Macauley(nextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + thisFabric);
            dPStrain2   = nextDGamma * contravariant2Covariant(R);
        end 
    end
    aCep2 = getElastoPlasticTangent(parms, thisSigma, nextDGamma, G, K, B, C, D, h, n, b);
    
    % update
    dSigma  =  (dSigma1 + dSigma2) / 2.0;
    dAlpha  =  (dAlpha1 + dAlpha2) / 2.0;
    dFabric =  (dFabric1 + dFabric2) / 2.0;
    dPStrain = (dPStrain1 + dPStrain2) / 2.0;
    nStress = nextStress + dSigma;
    nAlpha  = nextAlpha  + dAlpha;
    nFabric = nextFabric  + dFabric;
    
    p = one3 * getTrace(nextStress) + m_Presidual;
    if p < m_Presidual
        if dT == dT_min
            return;
        end
        dT = max(0.1 * dT, dT_min);
        continue;
    end
    stressNorm = getNormalContr(nextStress);
    % alphaNorm = getNormalContr(nextAlpha);
    temp0 = dSigma2 - dSigma1;
    if stressNorm < 0.5
        curStepError = getNormalContr(temp0);
    else
        curStepError = getNormalContr(temp0) / (2 * stressNorm);
    end
    fprintf('iteration: %d, error: %.5e, T: %.5e, dT: %.5e \n',iteration, curStepError, T, dT);
    if curStepError > TolE
        fprintf('--- unsuccessful increment\n ');
        q = max(0.8 * sqrt(TolE / curStepError), 0.1);
        if dT==dT_min
            nextElasticStrain = nextElasticStrain - dPStrain;
            nextStress = nStress;
            eta = sqrt(13.5) * getNormalContr(getDevPart(nextStress)) / getTrace(nextStress);
            nextAlpha  = curAlpha + 3.0 * (getDevPart(nextStress)/getTrace(nextStress) - getDevPart(curStress)/getTrace(curStress));
            T = T + dT;
        end
        dT = max(q * dT, dT_min);
    else
        fprintf('+++ successful increment\n ');
        nextElasticStrain = nextElasticStrain - dPStrain;
        nextStress = nStress;
        nextAlpha  = nAlpha;
        nextFabric = nFabric;
        [nextStress,nextAlpha,nextElasticStrain,aCep,aCep_Consistent] = stressCorrection(parms, curStress, curAlpha, curElasticStrain,...
        nextStress, nextAlpha, nextFabric,nextElasticStrain, nextVoidRatio, alpha_in, G,K,aCep, aCep_Consistent);
       
        T = T + dT;
        aCep_thisStep =   (aCep1 + aCep2) / 2.0;
        aCep_Consistent = aCep_thisStep * (aD * aCep_Consistent + T * mIImix);

        q = max(0.8 * sqrt(TolE / curStepError), 0.5);
        dT = max(q * dT, dT_min);
        dT = min(dT, 1 - T);
    end
    iteration = iteration + 1;
end
fn = getF(parms,nextStress,nextAlpha) % 1.938544
end


%%
function [nStress,nAlpha,nElasticStrain,naCep,naCep_Consistent] = stressCorrection(parms, curStress, curAlpha, curElasticStrain,...
          nextStress, nextAlpha, nextFabric,nextElasticStrain, nextVoidRatio, alpha_in, G, K, aCep, aCep_Consistent)
% Purpose: stressCorrection correct the stress state to the yield surface
% Arguments:
%                           I/O   Type
%  parms                     I    R    : model constants given by the user
%  curStress                 I    6x1  : current stress tensor, contravariant form
%  curAlpha                  I    6x1  : current back stress tensor, contravariant form
%  curElasticStrain          I    6x1  : current elastic strain tensor, contravariant form
%  nextStress                I    6x1  : next stress tensor, contravariant form
%  nextAlpha                 I    6x1  : next back stress tensor, contravariant form
%  nextFabric                I    6x1  : next fabric tensor, contravariant form
%  nextElasticStrain         I    6x1  : next elastic strain tensor, contravariant form
%  nextVoidRatio             I    R    : next void ratio
%  alpha_in                  I    R    : initial fabric anisotropy
% Local variables
one3 = 1.0/3.0;
two3 = 2.0/3.0;
mTolF = parms(18);
m_Pmin = parms(22);
m_Presidual = parms(23);
maxIter = 50;
mI1 = [1.0; 1.0; 1.0; 0.0; 0.0; 0.0];
p = one3 * getTrace(curStress)+m_Presidual;
if p < m_Pmin + m_Presidual
    p = m_Pmin + m_Presidual;
    nStress = p * mI1;
    nAlpha = nextAlpha;
    nElasticStrain = nextElasticStrain;
    naCep = aCep;
    naCep_Consistent = aCep_Consistent;
    return;
else
    % See if NextStress is outside yield surface
    fr = getF(parms, nextStress, nextAlpha);
    if abs(fr) < mTolF  % NextStress is inside yield surface
        nStress = nextStress;
        nAlpha = nextAlpha;
        nElasticStrain = nextElasticStrain;
        naCep = aCep;
        naCep_Consistent = aCep_Consistent;
        return;
    else
        nStress = nextStress;
        nAlpha = nextAlpha;
        for i = 1:maxIter
            fprintf('Stress state outside yield surface. Correction step = %d, yield function = %.5e \n', i, fr );
            devStress = getDevPart(nStress);
            aC = getElasticStiffness(G, K);
            aD = getCompliance(G, K);
            [n, d, b, cos3Theta, h, psi, alphaBtheta, alphaDtheta,b0, A_d, D, B, C, R] = getStateDependent(...
            parms, nStress, nAlpha, nextFabric, nextVoidRatio, alpha_in);
            dSigmaP = doubleDot4_2(aC, covariant2Contraviant(R));
            aBar = two3 * h * b;
            r = devStress / p ;
            dfrOverdSigma = n - one3 * doubleDot2_2_Contr(n, r) * mI1;
            dfrOverdAlpha = - p * n;
            lambda = fr / (doubleDot2_2_Contr(dfrOverdSigma, dSigmaP)-doubleDot2_2_Contr(dfrOverdAlpha, aBar));
            if  abs(getF(parms, nStress - lambda * dSigmaP, nAlpha + lambda * aBar)) < abs(fr)
                nStress = nStress - lambda * dSigmaP;
                nAlpha = nAlpha + lambda * aBar;
            else
                % lambda = fr / doubleDot2_2_Contr(dfrOverdSigma, dfrOverdSigma)
                % nStress = nStress - lambda * dfrOverdSigma;
                if abs(getF(parms, nStress - lambda * dfrOverdSigma, nAlpha)) < abs(fr)
                    nStress = nStress - lambda * dfrOverdSigma;
                else
                    disp('Warning! Stress correction failed.');
                    nElasticStrain = nextElasticStrain;
                    naCep = aCep;
                    naCep_Consistent = aCep_Consistent;
                    return;
                end
            end
            fr = getF(parms, nStress, nAlpha);
            if abs(fr) < mTolF
                nextStress = nStress;
                nextAlpha  = nAlpha;
                break;
            end
            if i == maxIter
                fprintf('Stile outside yield surface, use another method. Correction step = %d, yield function = %.5e \n', i, fr);
                if (getF(parms, curStress, nextAlpha) < mTolF)
                    dSigma = nextStress - curStress;
                    alpha_up  = 1.0;
                    alpha_mid = 0.5;
                    alpha_low = 0.0;
                    fr_old = getF(parms, curStress+alpha_mid * dSigma, nextAlpha);
                    for j = 1 : maxIter
                        if fr_old  < 0.0
                            alpha_down = alpha_mid;
                            alpha_mid = 0.5 * (alpha_up + alpha_down);
                        else
                            alpha_up = alpha_mid;
                            alpha_mid = 0.5 * (alpha_up + alpha_low);
                        end
                        fr_old = getF(parms, curStress+alpha_mid * dSigma, nextAlpha);
                        
                        if abs(fr_old) < mTolF
                            nStress = curStress + alpha_mid * dSigma;
                            break;
                        end
                        if jj == maxIter
                            disp('Warning! Stress correction failed by second method.');
                            break;
                        end
                    end
                else
                    nextStress = curStress;
                    nextAlpha  = curAlpha;
                end
            end
            p = one3 * getTrace(nextStress) + m_Presidual;
        end
    end
    nElasticStrain =  curElasticStrain + doubleDot4_2(aD, nextStress - curStress);
    naCep = getElastoPlasticTangent(parms, nextStress, nextAlpha, G, K, B, C, D, h, n, b);
    naCep_Consistent = naCep; % n means new, output of this function
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   SYMMETRIC TENSOR OPERATIONS  %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In all the functions below, contravariant means a stress-like tensor
% and covariant means a strain-like tensor


%%
function trace = getTrace(vector)
% Purpose: compute the trace of the input argument
if (length(vector) ~= 6)
    error('Error! getTrace requires vector of size(6)!');
else
    trace = sum(vector(1:3));
end
end


%%
function devpart = getDevPart(vector)
% Purpose: computes the deviatoric part of the input tensor
if (length(vector) ~= 6)
    error('Error! GetDevPart requires vector of size(6)!');
else
    one3 = 1.0/3.0;
    p = getTrace(vector);
    devpart = vector;
    devpart(1) = devpart(1) - one3 * p;
    devpart(2) = devpart(2) - one3 * p;
    devpart(3) = devpart(3) - one3 * p;
end
end


%%
function singleDot = singleDot(v1, v2)
% Purpose: computes v1.v2, v1 and v2 should be both in their "contravariant" form
if (length(v1) ~= 6) || (length(v2) ~= 6)
    error('Error! singleDot requires vector of size(6)!');
else
    singleDot = zeros(6,1);
    singleDot(1) = v1(1)*v2(1) + v1(4)*v2(4) + v1(6)*v2(6);
    singleDot(2) = v1(4)*v2(4) + v1(2)*v2(2) + v1(5)*v2(5);
    singleDot(3) = v1(6)*v2(6) + v1(5)*v2(5) + v1(3)*v2(3);
    singleDot(4) = 0.5*(v1(1)*v2(4) + v1(4)*v2(1) + v1(4)*v2(2) + v1(2)*v2(4) + v1(6)*v2(5) + v1(5)*v2(6));
    singleDot(5) = 0.5*(v1(4)*v2(6) + v1(6)*v2(4) + v1(2)*v2(5) + v1(5)*v2(2) + v1(5)*v2(3) + v1(3)*v2(5));
    singleDot(6) = 0.5*(v1(1)*v2(6) + v1(6)*v2(1) + v1(4)*v2(5) + v1(5)*v2(4) + v1(6)*v2(3) + v1(3)*v2(6));
end
end


%%
function result = doubleDot2_2_Contr(v1, v2)
% Purpose: computes v1:v2, v1 and v2 should be both in their "contravariant" form
if (length(v1) ~= 6) || (length(v2) ~= 6)
    error('Error! doubleDot2_2_Contr requires vector of size(6)!');
else
    result = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3) + 2.0*v1(4)*v2(4) + 2.0*v1(5)*v2(5) + 2.0*v1(6)*v2(6);
end
end


%%
function result = doubleDot2_2_Cov(v1, v2)
% Purpose: computes v1:v2, v1 and v2 should be both in their "covariant" form
if (length(v1) ~= 6) || (length(v2) ~= 6)
    error('Error! DoubleDot2_2_Cov requires vector of size(6)!');
else
    result = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3) - 0.5*v1(4)*v2(4) - 0.5*v1(5)*v2(5) - 0.5*v1(6)*v2(6);
end
end


%%
function result = doubleDot2_2_Mixed(v1, v2)
% Purpose: computes doubledot product for vector-vector arguments, one "covariant" and the other "contravariant"
if (length(v1) ~= 6) || (length(v2) ~= 6)
    error('Error! DoubleDot2_2_Mixed requires vector of size(6)!');
else
    result = 0.0;
    for i = 1:6
        result = result + v1(i)*v2(i);
    end
end
end


%%
function contravariant = toContravariant(tensor)
% Purpose: converts a tensor to its contravariant form, i.e. stress-like
if isequal(size(tensor), [3,3])
    % contravariant = zeros(6,1);
    contravariant = [tensor(1,1); tensor(2,2); tensor(3,3); tensor(1,2); tensor(2,3); tensor(1,3)];
else
    error('Error! ToContravariant requires 3x3 matrix!');
end
end


%%
function covariant = toCovariant(tensor)
% Purpose: converts a tensor to its covariant form, i.e. strain-like
if isequal(size(tensor), [3,3])
    covariant = [tensor(1,1); tensor(2,2); tensor(3,3); 2.0*tensor(1,2); 2.0*tensor(2,3); 2.0*tensor(1,3)];
else
    error('Error! ToCovariant requires 3x3 matrix!');
end
end


%%
function covariant = contravariant2Covariant(v1)
% Purpose: converts a vector to its covariant form, i.e. stress-like -> strain-like
if (length(v1) ~= 6) 
    error('Error! Contravariant2Covariant requires vector of size(6)!');
else
    covariant = v1;
    covariant(4)=v1(4) * 2.0; 
    covariant(5)=v1(5) * 2.0;
    covariant(6)=v1(6) * 2.0;
end
end


%%
function contravariant = covariant2Contraviant (v1)
% Purpose: converts a vector to its contravariant form, i.e. strain-like -> stress-like
if (length(v1) ~= 6) 
    error('Error! covariant2contravariant requires vector of size(6)!');
else
    one2 = 0.5;
    contravariant = v1;
    contravariant(4)=v1(4) * one2; 
    contravariant(5)=v1(5) * one2;
    contravariant(6)=v1(6) * one2;
end
end


%%
function result = doubleDot4_2(m1, v1)
% Purpose: computes doubledot product for matrix-vector arguments
if (length(v1) ~= 6)
    error('Error! doubleDot4_2 requires vector of size(6)!');
elseif (size(m1,1) ~= 6) || (size(m1,2) ~= 6)
    error('Error! doubleDot4_2 requires 6x6 matrix!');
else
    result = m1*v1;
end
end


%%
function result = MacauleyIndex(x)
% Purpose: computes Macauley index
if x >= 0
    result = 1;
else
    result = 0;
end
end


%%
function result = Macauley(x)
% Purpose: computes Macauley function
if x >= 0
    result = x;
else
    result = 0;
end
end


%%
function result = getNormalContr(vector)
% Purpose: compute contravariant (stress-like) norm of input 6x1 tensor

if (length(vector) ~= 6)
    error('Error! getNormalContr requires vector of size(6)!');
else
    result = sqrt(doubleDot2_2_Contr(vector,vector));
end
end


%%
function [p,q] = inv_s(sigma)
% Purpose: Invariants of the stress tensor sigma 
%        p=tr(sig)/3
%        q=sqrt((3/2)*tr(s*s))
fact=1/3;
m=[1;  1;  1;  0;  0;  0];
M=[1,  0,  0,  0,  0,  0
    0,  1,  0,  0,  0,  0
    0,  0,  1,  0,  0,  0
    0,  0,  0,  2,  0,  0
    0,  0,  0,  0,  2,  0
    0,  0,  0,  0,  0,  2];

p=fact*m'*sigma;

s=sigma-p*m;
sn=s'*(M*s);
if sigma(3)>sigma(1)
    q=sqrt(3*sn/2); % when z axis is the major principal axis
else
    q=-sqrt(3*sn/2);
end
end


%%
function [ev,eq] = inv_e(eps)
% Purpose: Invariants of the strain tensor eps
%        ev=tr(eps)
%        q=sqrt((2/3)*tr(e*e))
%        J2e=e_ij*e_ji a
fact=1/3;
m=[1;  1;  1;  0;  0;  0];
Minv=[1,  0,  0,   0,   0,   0;
        0,  1,  0,   0,   0,   0;
        0,  0,  1,   0,   0,   0;
        0,  0,  0, 0.5,   0,   0;
        0,  0,  0,   0, 0.5,   0;
        0,  0,  0,   0,   0, 0.5];
ev=m'*eps;
e=eps-(ev*fact)*m;
J2e=e'*(Minv*e);
if eps(3)>eps(1)
    eq=sqrt(2*J2e/3); % when z axis is the major principal axis
else
    eq=-sqrt(2*J2e/3);
end
end


%%
function draw_picture(parms,output)
% Purpose: draw the figures
% Local variables
one2 = 0.5;
stress = output(:,7:12);
fabric = output(:,13:18);
strain = output(:,25:30);
e = output(:,31);
yf_update_flag = output(:,32);
cycle = output(:,34);
%  getTrace(stress(100,:))
for i = 1:1:size(stress,1)
    % p(i) = sum(stress(i,:))/3.0;
    [p(i),q(i)] = inv_s(stress(i,:)');
    [epsv(i),epsq(i)] = inv_e(strain(i,:)');
end
figure(1)

% q  VS epsq, axis direction: x_3
subplot(2,3,1)
plot(epsq,q,'b-')
xlabel('deviatoric strain, \epsilon_q (%)')
ylabel('deviator stress, q (kPa)')
grid on
hold on

% q VS p 
subplot(2,3,2)
plot(p,q,'b-')
xlabel('mean effective stress, p (kPa)')
ylabel('deviator stress, q (kPa)')
grid on
hold on
% axis equal

% epsv VS epsq
subplot(2,3,3)
plot(epsq,epsv,'b-')
xlabel('deviatoric strain, \epsilon_q (%)')
ylabel('volumetric strain, \epsilon_v (%)')
grid on
hold on

% p VS epsv
subplot(2,3,4)
plot(p,epsv,'b-')
xlabel('mean effective stress, p (kPa)')
ylabel('volumetric strain, \epsilon_v (%)')
grid on
hold on

% e VS p
subplot(2,3,5)
plot(p,e,'b-')
xlabel('mean effective stress, p (kPa)')
ylabel('void ratio e')
grid on
hold on

figure(2)
clf
% q vs p and yield surface when monotonic loading
plot(p,q,'b-')
xlabel('mean effective stress, p (kPa)')
ylabel('deviator stress, q (kPa)')
grid on
hold on
x_plot = 0:1/100:1;
m_m = parms(8);
p_max_plot = max(p()) * 1.5;
p_plot = x_plot * p_max_plot;
yf_update_flag_plot = yf_update_flag(1:end);
% plot yield surface during the loading
number = min(floor(size(p,2)/3), 30);
number = min(number, size(p,2));
points_plot = round(linspace(2, size(p,2),number));
points_plot = unique([2,points_plot,size(p,2)], 'stable'); % add the second and last points
for i = 1:1:size(points_plot,2)
    if yf_update_flag_plot(points_plot(i)) == 0
        
        eta_upper = q(points_plot(i))/p(points_plot(i));
        eta_lower = eta_upper - 2 * m_m;
        q_upper = eta_upper * p_plot;
        q_lower = eta_lower * p_plot;
        plot(p_plot,q_upper,'r-')
        plot(p_plot,q_lower,'r-')
        plot(p(points_plot(i)),q(points_plot(i)),'ro')
    end
end
% plot initial yield surface 
q_upper = m_m * p_plot;
q_lower = -m_m * p_plot;
plot(p_plot,q_upper,'k-')
plot(p_plot,q_lower,'k-')
plot(p(1),q(1),'ko')
end
    