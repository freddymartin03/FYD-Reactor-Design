clear 

%Initial Conditions 
nMo = 1036172;
nTo = 1036172;
nAo = 5206.894;
nPXo = 0;
nOXo = 0;
nMXo = 0;
nBo = 21146.36;
nEo = 0;
nWo = 0;
Ex1o = 0;
Ex2o = 0;
Ex3o = 0;
Ex4o = 0;
PXtar = 1036172;

%revenues/costs [$/kmol]
cT=62.6552/0.98;
cM=12.17596/0.995;
cPX=139.0696;
cOX=71.1272;
cMX=74.312;
cB=64.83462;
cCat=1.1*0.55*1310000;
cE=1411.2/1055*3.8;

%number of iterations ODE solver does
opts = odeset('MaxStep',1e-2);

%Time for which integration is done over
mspan=[0 500];

To = 820;
Po = 6;
% Initial conditions vector
xo=[nTo, nMo, nPXo, nOXo, nMXo, nBo, nEo, nAo, nWo, Ex1o, Ex2o, Ex3o, Ex4o, To, Po];
%ODE solver
[v,x]=ode15s(@Reactor,mspan ,xo,opts);

% defining the output of each column of x from ODE solver
    nT = x(:,1);
    nM = x(:,2);
    nPX = x(:,3);
    nOX = x(:,4);
    nMX = x(:,5);
    nB = x(:,6);
    nE = x(:,7);
    nA = x(:,8);
    nW = x(:,9);
    Ex1 = x(:,10);
    Ex2 = x(:,11);
    Ex3 = x(:,12);
    Ex4 = x(:,13);
    frac = PXtar./nPX;
    Efactor = (nOX+nMX+nB+nE)./nPX;
    nTotal = nT+nM+nPX+nOX+nMX+nB+nE+nA+nW;
        
Conv = (nTo-nT)./nTo;
EP2a = PXtar*cPX+(nOX.*cOX+nMX.*cMX+nB.*cB+nE.*cE-(nTo-nT).*cT-(nMo-nM).*cM).*PXtar./nPX;
EP2b = PXtar*cPX+(nOX.*cOX+nMX.*cMX+nB.*cB+nE.*cE-nTo.*cT-nMo.*cM).*PXtar./nPX;
EP2c = PXtar*cPX+(nOX.*cOX+nMX.*cMX+nB.*cB+nE.*cE-nTo.*cT-(nMo-nM).*cM).*PXtar./nPX;
EP2d = PXtar*cPX+(nOX.*cOX+nMX.*cMX+nB.*cB+nE.*cE-(nTo-nT).*cT-nMo.*cM).*PXtar./nPX;
EP2anb = PXtar*cPX+(nOX.*cOX+nMX.*cMX+nE.*cE-(nTo-nT).*cT-(nMo-nM).*cM).*PXtar./nPX;
EP2bnb = PXtar*cPX+(nOX.*cOX+nMX.*cMX+nE.*cE-nTo.*cT-nMo.*cM).*PXtar./nPX;
EP2cnb = PXtar*cPX+(nOX.*cOX+nMX.*cMX+nE.*cE-nTo.*cT-(nMo-nM).*cM).*PXtar./nPX;
EP2dnb = PXtar*cPX+(nOX.*cOX+nMX.*cMX+nE.*cE-(nTo-nT).*cT-nMo.*cM).*PXtar./nPX;

figure(5)
hold on 
plot(Conv, nT./nTotal, 'LineWidth', 1.4)
plot(Conv, nM./nTotal, 'LineWidth', 1.4)
plot(Conv, nPX./nTotal, 'LineWidth', 1.4)
plot(Conv, nOX./nTotal, 'LineWidth', 1.4)
plot(Conv, nMX./nTotal, 'LineWidth', 1.4)
plot(Conv, nB./nTotal, 'LineWidth', 1.4)
plot(Conv, nE./nTotal, 'LineWidth', 1.4)
plot(Conv, nA./nTotal, 'LineWidth', 1.4, 'Color', [1 0 1])
plot(Conv, nW./nTotal, 'LineWidth', 1.4, 'Color', [0 0 1])
xlim([0.1 0.6])
xlabel('Conversion')
ylabel('Mole Fraction')
lgd = legend('Toluene', 'Methanol', 'p-Xylene', 'o-Xylene', 'm-Xylene', 'Benzene', 'Ethylene', 'Acetonitrile', 'Water');
lgd.NumColumns = 2;
title(lgd, 'Component')
hold off

To = 730;
Po = 1;
nTo = 1036172*4;
nBo = 21146.36*4;
% Initial conditions vector
xo=[nTo, nMo, nPXo, nOXo, nMXo, nBo, nEo, nAo, nWo, Ex1o, Ex2o, Ex3o, Ex4o, To, Po];
%ODE solver
[v,x]=ode15s(@Reactor,mspan ,xo,opts);

% defining the output of each column of x from ODE solver
    nT = x(:,1);
    nM = x(:,2);
    nPX = x(:,3);
    nOX = x(:,4);
    nMX = x(:,5);
    nB = x(:,6);
    nE = x(:,7);
    nA = x(:,8);
    nW = x(:,9);
    Ex1 = x(:,10);
    Ex2 = x(:,11);
    Ex3 = x(:,12);
    Ex4 = x(:,13);
    frac4 = PXtar./nPX;
    Efactor4=(nOX+nMX+nB+nE)./nPX;


Conv4 = (nTo-nT)./nTo;
EP2aalt = PXtar*cPX+(nOX.*cOX+nMX.*cMX+nB.*cB+nE.*cE-(nTo-nT).*cT-(nMo-nM).*cM).*PXtar./nPX;
EP2balt = PXtar*cPX+(nOX.*cOX+nMX.*cMX+nB.*cB+nE.*cE-nTo.*cT-nMo.*cM).*PXtar./nPX;
EP2calt = PXtar*cPX+(nOX.*cOX+nMX.*cMX+nB.*cB+nE.*cE-nTo.*cT-(nMo-nM).*cM).*PXtar./nPX;
EP2dalt = PXtar*cPX+(nOX.*cOX+nMX.*cMX+nB.*cB+nE.*cE-(nTo-nT).*cT-nMo.*cM).*PXtar./nPX;
EP2aaltnb = PXtar*cPX+(nOX.*cOX+nMX.*cMX+nE.*cE-(nTo-nT).*cT-(nMo-nM).*cM).*PXtar./nPX;
EP2baltnb = PXtar*cPX+(nOX.*cOX+nMX.*cMX+nE.*cE-nTo.*cT-nMo.*cM).*PXtar./nPX;
EP2caltnb = PXtar*cPX+(nOX.*cOX+nMX.*cMX+nE.*cE-nTo.*cT-(nMo-nM).*cM).*PXtar./nPX;
EP2daltnb = PXtar*cPX+(nOX.*cOX+nMX.*cMX+nE.*cE-(nTo-nT).*cT-nMo.*cM).*PXtar./nPX;

figure(1)
hold on
plot(Conv,EP2a/1000000)
plot(Conv,EP2b/1000000)
plot(Conv,EP2c/1000000)
plot(Conv,EP2d/1000000)
xlim([0.1  0.6])
ylim([-100 100])
xlabel('Conversion')
ylabel('Economic Potential (Million $/year')
legend('Base','Alt1','Alt2','Alt3')
hold off

figure(2)
hold on
plot(Conv,EP2anb/1000000)
plot(Conv,EP2bnb/1000000)
plot(Conv,EP2cnb/1000000)
plot(Conv,EP2dnb/1000000)
xlim([0.1  0.6])
ylim([-100 100])
xlabel('Conversion')
ylabel('Economic Potential (Million $/year')
legend('Base','Alt1','Alt2','Alt3')
hold off

figure(3)
hold on
plot(Conv4,EP2aalt/1000000)
plot(Conv4,EP2balt/1000000)
plot(Conv4,EP2calt/1000000)
plot(Conv4,EP2dalt/1000000)
xlim([0.1  0.6])
ylim([-200 500])
xlabel('Conversion')
ylabel('Economic Potential (Million $/year')
legend('Base','Alt1','Alt2','Alt3')
hold off

figure(4)
hold on
plot(Conv4,EP2aaltnb/1000000)
plot(Conv4,EP2baltnb/1000000)
plot(Conv4,EP2caltnb/1000000)
plot(Conv4,EP2daltnb/1000000)
xlim([0.1  0.6])
xlabel('Conversion')
ylabel('Economic Potential (Million $/year')
legend('Base','Alt1','Alt2','Alt3')
hold off



function dxdv = Reactor(v,x)

nT = x(1);
nM = x(2);
nPX = x(3);
nOX = x(4);
nMX = x(5);
nB = x(6);
nE = x(7);
nA = x(8);
nW = x(9);
Ex1 = x(10);
Ex2 = x(11);
Ex3 = x(12);
Ex4 = x(13);
T = x(14);
P = x(15);

% General Data
Catden = 1310000;
void = 0.45;
k1_0=464*8;
k2_0=249*8;
k3_0=4.4*8;
k4_0 = 163*8;
Ea1 = 63.5*1000;
Ea2 = 68.5*1000;
Ea3 = 23.6*1000;
Ea4 = 51.2*1000;
n1 = 0.728;
n2 = 1;
n3 = 1.55;
n4 = 1;
R = 8.3145;

%rate constant
k1 = k1_0*exp(-Ea1/(R*T));
k2 = k2_0*exp(-Ea2/(R*T));
k3 = k3_0*exp(-Ea3/(R*T));
k4 = k4_0*exp(-Ea4/(R*T));

%overall material balance
nTotal = nT + nM + nPX + nOX + nMX + nB + nE + nA + nW;

% Partial Pressures
PT = nT .* P ./ nTotal;
PM = nM.* P ./ nTotal;
PPX = nPX.* P ./ nTotal;

% Rate Expressions
r1 = n1*k1*PT*PM;
r2 = n2*k2*PT;
r3 = n3*k3*PPX;
r4 = n4*k4*PPX;

% list of ODEs
%material balances
dnTdv = (r3-r2-r1)*(1-void)*Catden;
dnMdv = -r1*(1-void)*Catden;
dnPXdv = (r1+0.5*r2-r3-r4)*(1-void)*Catden;
dnOXdv = (0.5*r4)*(1-void)*Catden;
dnMXdv = (0.5*r4)*(1-void)*Catden;
dnBdv = (0.5*r2)*(1-void)*Catden;
dnEdv = (0.5*r3)*(1-void)*Catden;
dnAdv = 0;
dnWdv = (0.5*r1)*(1-void)*Catden;

%extent
dEx1dv= r1*(1-void)*Catden;
dEx2dv= r2*(1-void)*Catden;
dEx3dv= r3*(1-void)*Catden;
dEx4dv= r4*(1-void)*Catden;

%T and P constant
dTdv = 0;
dPdv = 0;

%vector of differentials for ODE solver
dxdv=[dnTdv;dnMdv;dnPXdv;dnOXdv;dnMXdv;dnBdv;dnEdv;dnAdv;dnWdv;dEx1dv;dEx2dv;dEx3dv;dEx4dv;dTdv;dPdv];
end

