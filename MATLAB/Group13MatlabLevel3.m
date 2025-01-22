clear 
clc

%Initial Conditions 
nMo = 1036172; 
nAo = 5206.8940;
nPXo = 0;
nOXo = 0;
nMXo = 0;
nEo = 0;
nWo = 0;
Ex1o = 0;
Ex2o = 0;
Ex3o = 0;
Ex4o = 0;
PXtar = 1036172;

%revenues/costs (dollars/kg) 
cT=62.6552/0.98;
cM=12.17596/0.995;
cPX=139.0696;
cOX=71.1272;
cMX=74.312;
cB=64.83462;
cCat= 1.1; 
cfuel = 3.8; % $ per 10^6 btu
cE=1411.2/1055*3.8;

%Cp values for feed preheat
fcpT = 231.4185;
fcpM = 79.86975;
fcpPX = 275.886;
fcpOX = 277.5705;
fcpMX = 275.931;
fcpB = 187.728;
fcpE = 83.5651;
fcpA = 89.3958;
fcpW = 39.59435;

%Hvap values for feed preheat
HT = 37300;
HM = 37600;
HPX = 42000;
HOX = 42000;
HMX = 42000;
HB = 33900;
HE = 13600;
HA = 33400;
HW = 40600;

%Indices and other constants 
MS = 1800; 
Fc = 1; 
voidage = 0.45; 
catdensity = 1310000/1000; % kg/m^3 
k = 10;

opts = odeset('MaxStep',1e-2);

%Time for which integration is done over
mspan=[0 2];
Po = 6;
nTo = 1036172;
nBo = 21146.36;
mat=zeros(900,3);
for b=1:900
    To=729.9+0.1*b; %range of values for bed voidage

    xo=[nTo, nMo, nPXo, nOXo, nMXo, nBo, nEo, nAo, nWo, Ex1o, Ex2o, Ex3o, Ex4o, To, Po];
    [v,x]=ode15s(@AdiabaticReactor,mspan ,xo,opts);
    
    T = x(:,14);
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
    Conv = (nTo-nT)./nTo;

   
    EP2 = PXtar.*cPX+(nOX*cOX+nMX*cMX+nB*cB+nE*cE-(nTo-nT)*cT-(nMo-nM)*cM).*PXtar./nPX;
    %+nB*cB
    %EP3 calculations
    %reactor cost
    volume = v*PXtar./nPX;
    D = ((volume.*(4/(3.1416*k))).^(1/3))*3.28084; % has to be in ft.
    L = k.*D; 
    BMC = (MS/280) .*101.9.* (D.^(1.066)).*(L.^(0.802)).*(2.18 + Fc);
    Catcost = volume.*(1-voidage).*catdensity.*cCat ;

    %furnace cost
    Qheat = ((nTo*fcpT+nMo*fcpM+nBo*fcpB+nAo*fcpA)*(To-293)+nTo*HT+nMo*HM+nBo*HB+nAo*HA)*PXtar./nPX;
    Qheatbtu = ((Qheat/1.055)/10^6)/8000;
    fuelcost = Qheatbtu*cfuel*8000;
    furncost = (MS/280)*(5520)*Qheatbtu.^(0.85)*(1.27+Fc);

    TCI = 2.09 .*(BMC +furncost)+ Catcost;
    TOC = fuelcost;
    
    TC = TCI./3 + TOC;
    
    EP3 = EP2 - TC ;
    [T uidx] = unique(T);
    Conv=Conv(uidx);
    EP3 = EP3(uidx);
    Qheat = Qheat(uidx);
    mat(b,1)=interp1(T,Conv,820);
    mat(b,2)=interp1(T,EP3,820);
    mat(b,3)=interp1(T,Qheat,820);
end


%Isothermal reactor
    H1 = -69825.5;
    H2 =325.5;
    H3 = 58115.425;
    H4 = 1017.5;

nTo = 1036172;
nBo = 21146.36;
To=820;
Po = 6;
xo=[nTo, nMo, nPXo, nOXo, nMXo, nBo, nEo, nAo, nWo, Ex1o, Ex2o, Ex3o, Ex4o, To, Po];
[v,x]=ode15s(@IsothermalReactor,mspan ,xo,opts);

T = x(:,14);
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
nTotal = nT+nPX+nOX+nMX+nB+nE+nA+nW;
fracI = nPX./nTotal;
Convbhs = (nTo-nT)./nTo;

EP2 = PXtar.*cPX+(nOX*cOX+nMX*cMX+nB*cB+nE*cE-(nTo-nT)*cT-(nMo-nM)*cM).*PXtar./nPX;

QR=(H1*Ex1+H2*Ex2+H3*Ex3+H4*Ex4)*PXtar./nPX;

volume = v*PXtar./nPX;
D = ((volume.*(4/(3.1416*k))).^(1/3))*3.28084; % has to be in ft.
L = k.*D; 
BMC = (MS/280) .*101.9.* (D.^(1.066)).*(L.^(0.802)).*(2.18 + Fc);
Catcost = volume.*(1-voidage).*catdensity.*cCat ;

%heat exchanger cost
deltaT=((To-308)-(To-288))/(log((To-308)/(To-288)));
Area = -QR*(1000/(8000*3600))/(deltaT*114);
QBTU = ((-QR/1.055)/10^6)/8000;
A = Area*10.764;
HXC=(MS/280)*101.3*A.^(0.65)*(2.29+0.85);
water = -QR/(1000*4.2*20);
watcost=8.12*water/1000;

%furnace cost
Qheat = ((nTo*fcpT+nMo*fcpM+nBo*fcpB+nAo*fcpA)*(To-293)+nTo*HT+nMo*HM+nBo*HB+nAo*HA)*PXtar./nPX;
Qheatbtu = ((Qheat/1.055)/10^6)/8000;
fuelcost = Qheatbtu*cfuel*8000;
furncost = (MS/280)*(5520)*Qheatbtu.^(0.85)*(1.27+Fc);

TCI = 2.09 .*(BMC+HXC+furncost) + Catcost;

TOC = watcost+fuelcost;

TC = TCI./3 + TOC;

EP3bhs = EP2 - TC ;

% [ConvIun uidx] = unique(ConvIun);
% fracI = fracI(uidx);
% molefracparaI = interp1(ConvIun,fracI,0.2)


%no benzene high selectivity cases
matnb=zeros(900,2);
for b=1:900
    To=729.9+0.1*b; %range of values for bed voidage

    xo=[nTo, nMo, nPXo, nOXo, nMXo, nBo, nEo, nAo, nWo, Ex1o, Ex2o, Ex3o, Ex4o, To, Po];
    [v,x]=ode15s(@AdiabaticReactor,mspan ,xo,opts);
    
    T = x(:,14);
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
    Conv = (nTo-nT)./nTo;

   
    EP2 = PXtar.*cPX+(nOX*cOX+nMX*cMX+nE*cE-(nTo-nT)*cT-(nMo-nM)*cM).*PXtar./nPX;
    %+nB*cB
    %EP3 calculations
    %reactor cost
    volume = v*PXtar./nPX;
    D = ((volume.*(4/(3.1416*k))).^(1/3))*3.28084; % has to be in ft.
    L = k.*D; 
    BMC = (MS/280) .*101.9.* (D.^(1.066)).*(L.^(0.802)).*(2.18 + Fc);
    Catcost = volume.*(1-voidage).*catdensity.*cCat ;

    %furnace cost
    Qheat = ((nTo*fcpT+nMo*fcpM+nBo*fcpB+nAo*fcpA)*(To-293)+nTo*HT+nMo*HM+nBo*HB+nAo*HA)*PXtar./nPX;
    Qheatbtu = ((Qheat/1.055)/10^6)/8000;
    fuelcost = Qheatbtu*cfuel*8000;
    furncost = (MS/280)*(5520)*Qheatbtu.^(0.85)*(1.27+Fc);

    TCI = 2.09 .*(BMC +furncost)+ Catcost;
    TOC = fuelcost;
    
    TC = TCI./3 + TOC;
    
    EP3 = EP2 - TC ;
    [T uidx] = unique(T);
    Conv=Conv(uidx);
    EP3 = EP3(uidx);
    matnb(b,1)=interp1(T,Conv,820);
    matnb(b,2)=interp1(T,EP3,820);
end


%Isothermal reactor
    H1 = -69825.5;
    H2 =325.5;
    H3 = 58115.425;
    H4 = 1017.5;
    
nTo = 1036172;
nBo = 21146.36;
To=820;
Po = 6;
xo=[nTo, nMo, nPXo, nOXo, nMXo, nBo, nEo, nAo, nWo, Ex1o, Ex2o, Ex3o, Ex4o, To, Po];
[v,x]=ode15s(@IsothermalReactor,mspan ,xo,opts);

T = x(:,14);
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
nTotal = nT+nPX+nOX+nMX+nB+nE+nA+nW;
fracI = nPX./nTotal;
Convbhsnb = (nTo-nT)./nTo;


EP2 = PXtar.*cPX+(nOX*cOX+nMX*cMX+nE*cE-(nTo-nT)*cT-(nMo-nM)*cM).*PXtar./nPX;

QR=(H1*Ex1+H2*Ex2+H3*Ex3+H4*Ex4)*PXtar./nPX;

volume = v*PXtar./nPX;
D = ((volume.*(4/(3.1416*k))).^(1/3))*3.28084; % has to be in ft.
L = k.*D; 
BMC = (MS/280) .*101.9.* (D.^(1.066)).*(L.^(0.802)).*(2.18 + Fc);
Catcost = volume.*(1-voidage).*catdensity.*cCat ;

%heat exchanger cost
deltaT=((To-308)-(To-288))/(log((To-308)/(To-288)));
Area = -QR*(1000/(8000*3600))/(deltaT*114);
QBTU = ((-QR/1.055)/10^6)/8000;
A = Area*10.764;
HXC=(MS/280)*101.3*A.^(0.65)*(2.29+0.85);
water = -QR/(1000*4.2*20);
watcost=8.12*water/1000;

%furnace cost
Qheat = ((nTo*fcpT+nMo*fcpM+nBo*fcpB+nAo*fcpA)*(To-293)+nTo*HT+nMo*HM+nBo*HB+nAo*HA)*PXtar./nPX;
Qheatbtu = ((Qheat/1.055)/10^6)/8000;
fuelcost = Qheatbtu*cfuel*8000;
furncost = (MS/280)*(5520)*Qheatbtu.^(0.85)*(1.27+Fc);

TCI = 2.09 .*(BMC+HXC+furncost) + Catcost;

TOC = watcost+fuelcost;

TC = TCI./3 + TOC;

EP3bhsnb = EP2 - TC ;

%base adiabatic low selectivity
    opts = odeset('MaxStep',1e-1);

%Time for which integration is done over
mspan=[0 1000];
Po = 1;
nTo = 1036172*4;
nBo = 21146.36*4;
matls=zeros(90,3);
for b=1:90
    To=730+b; %range of values for bed voidage

    xo=[nTo, nMo, nPXo, nOXo, nMXo, nBo, nEo, nAo, nWo, Ex1o, Ex2o, Ex3o, Ex4o, To, Po];
    [v,x]=ode15s(@AdiabaticReactor,mspan ,xo,opts);
    
    T = x(:,14);
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
    Conv = (nTo-nT)./nTo;

   
    EP2 = PXtar.*cPX+(nOX*cOX+nMX*cMX+nB*cB+nE*cE-(nTo-nT)*cT-(nMo-nM)*cM).*PXtar./nPX;
    %+nB*cB
    %EP3 calculations
    %reactor cost
    volume = v*PXtar./nPX;
    D = ((volume.*(4/(3.1416*k))).^(1/3))*3.28084; % has to be in ft.
    L = k.*D; 
    BMC = (MS/280) .*101.9.* (D.^(1.066)).*(L.^(0.802)).*(2.18 + Fc);
    Catcost = volume.*(1-voidage).*catdensity.*cCat ;

    %furnace cost
    Qheat = ((nTo*fcpT+nMo*fcpM+nBo*fcpB+nAo*fcpA)*(To-293)+nTo*HT+nMo*HM+nBo*HB+nAo*HA)*PXtar./nPX;
    Qheatbtu = ((Qheat/1.055)/10^6)/8000;
    fuelcost = Qheatbtu*cfuel*8000;
    furncost = (MS/280)*(5520)*Qheatbtu.^(0.85)*(1.27+Fc);

    TCI = 2.09 .*(BMC +furncost)+ Catcost;
    TOC = fuelcost;
    
    TC = TCI./3 + TOC;
    
    EP3 = EP2 - TC ;
    [T uidx] = unique(T);
    Conv=Conv(uidx);
    EP3 = EP3(uidx);
    Qheat = Qheat(uidx);
    matls(b,1)=interp1(T,Conv,730);
    matls(b,2)=interp1(T,EP3,730);
    matls(b,3)=interp1(T,Qheat,820);
end

%%% alt3 isothermal high selec

mspan=[0 2];
%Isothermal reactor
    H1 = -69825.5;
    H2 =325.5;
    H3 = 58115.425;
    H4 = 1017.5;

    opts = odeset('MaxStep',1e-2);
nTo = 1036172;
nBo = 21146.36;
To=820;
Po = 6;
xo=[nTo, nMo, nPXo, nOXo, nMXo, nBo, nEo, nAo, nWo, Ex1o, Ex2o, Ex3o, Ex4o, To, Po];
[v,x]=ode15s(@IsothermalReactor,mspan ,xo,opts);

T = x(:,14);
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
nTotal = nT+nPX+nOX+nMX+nB+nE+nA+nW;
fracI = nPX./nTotal;
Convahs = (nTo-nT)./nTo;


EP2 = PXtar.*cPX+(nOX*cOX+nMX*cMX+nB*cB+nE*cE-(nTo-nT)*cT-nMo*cM).*PXtar./nPX;

QR=(H1*Ex1+H2*Ex2+H3*Ex3+H4*Ex4)*PXtar./nPX;

volume = v*PXtar./nPX;
D = ((volume.*(4/(3.1416*k))).^(1/3))*3.28084; % has to be in ft.
L = k.*D; 
BMC = (MS/280) .*101.9.* (D.^(1.066)).*(L.^(0.802)).*(2.18 + Fc);
Catcost = volume.*(1-voidage).*catdensity.*cCat ;

%heat exchanger cost
deltaT=((To-308)-(To-288))/(log((To-308)/(To-288)));
Area = -QR*(1000/(8000*3600))/(deltaT*114);
QBTU = ((-QR/1.055)/10^6)/8000;
A = Area*10.764;
HXC=(MS/280)*101.3*A.^(0.65)*(2.29+0.85);
water = -QR/(1000*4.2*20);
watcost=8.12*water/1000;

%furnace cost
Qheat = ((nTo*fcpT+nMo*fcpM+nBo*fcpB+nAo*fcpA)*(To-293)+nTo*HT+nMo*HM+nBo*HB+nAo*HA)*PXtar./nPX;
Qheatbtu = ((Qheat/1.055)/10^6)/8000;
fuelcost = Qheatbtu*cfuel*8000;
furncost = (MS/280)*(5520)*Qheatbtu.^(0.85)*(1.27+Fc);

TCI = 2.09 .*(BMC+HXC+furncost) + Catcost;

TOC = watcost+fuelcost;

TC = TCI./3 + TOC;

EP3ahs = EP2 - TC ;


%%%% alt3 adiabatic low selec
%base adiabatic low selectivity
    opts = odeset('MaxStep',1e-1);

%Time for which integration is done over
mspan=[0 1000];
Po = 1;
nTo = 1036172*4;
nBo = 21146.36*4;
matals=zeros(90,2);
for b=1:90
    To=730+b; %range of values for bed voidage

    xo=[nTo, nMo, nPXo, nOXo, nMXo, nBo, nEo, nAo, nWo, Ex1o, Ex2o, Ex3o, Ex4o, To, Po];
    [v,x]=ode15s(@AdiabaticReactor,mspan ,xo,opts);
    
    T = x(:,14);
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
    Conv = (nTo-nT)./nTo;

   
    EP2 = PXtar.*cPX+(nOX*cOX+nMX*cMX+nB*cB+nE*cE-(nTo-nT)*cT-nMo*cM).*PXtar./nPX;
    
    %EP3 calculations
    %reactor cost
    volume = v*PXtar./nPX;
    D = ((volume.*(4/(3.1416*k))).^(1/3))*3.28084; % has to be in ft.
    L = k.*D; 
    BMC = (MS/280) .*101.9.* (D.^(1.066)).*(L.^(0.802)).*(2.18 + Fc);
    Catcost = volume.*(1-voidage).*catdensity.*cCat ;

    %furnace cost
    Qheat = ((nTo*fcpT+nMo*fcpM+nBo*fcpB+nAo*fcpA)*(To-293)+nTo*HT+nMo*HM+nBo*HB+nAo*HA)*PXtar./nPX;
    Qheatbtu = ((Qheat/1.055)/10^6)/8000;
    fuelcost = Qheatbtu*cfuel*8000;
    furncost = (MS/280)*(5520)*Qheatbtu.^(0.85)*(1.27+Fc);

    TCI = 2.09 .*(BMC +furncost)+ Catcost;
    TOC = fuelcost;
    
    TC = TCI./3 + TOC;
    
    EP3 = EP2 - TC ;
    [T uidx] = unique(T);
    Conv=Conv(uidx);
    EP3 = EP3(uidx);
    matals(b,1)=interp1(T,Conv,730);
    matals(b,2)=interp1(T,EP3,730);
end

figure(1)
hold on
plot(mat(:,1),mat(:,2)/1000000)
plot(Convbhs,EP3bhs/1000000)
plot(matnb(:,1),matnb(:,2)/1000000)
plot(Convbhsnb,EP3bhsnb/1000000)
plot(matls(:,1),matls(:,2)/1000000)
plot(Convahs,EP3ahs/1000000)
plot(matals(:,1),matals(:,2)/1000000)
hold off
legend('base, adiabatic, high selec, sell Benzene','base, isothermal, high selec, sell Benzene','base, adiabatic, high selec, waste Benzene','base, isothermal, high selec, waste Benzene','base, adiabatic, low selec, sell Benzene','alt3, isothermal, high selec, sell Benzene','alt3, adiabatic, low selec, sell Benzene')
ylim([0 250])
xlim([0.1 0.6])
xlabel('Conversion')
ylabel('EP3 (million$/year)')


%% Function: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dxdv = AdiabaticReactor(v,x) 

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

    % %Heat Capacities Data (entered at 6 atm): to be entered: 
    CpT =230.5145;
    CpM =78.6352;
    CpPX =274.1495;
    CpOX =275.6325;
    CpMX=274.1535;
    CpB =187.458;
    CpE =83.01595;
    CpA =88.21465;
    CpW = 38.61675;

    % % Heat of Reaction Data:
    H1 = -69825.5;
    H2 =325.5;
    H3 = 58115.425;
    H4 = 1017.5;
        
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

    %T change
    QR = (r1.*H1 + r2.*H2 + r3.*H3 + r4.*H4)*(1-void)*Catden;
    Cpchange = nT.*CpT + nM.*CpM + nPX.*CpPX + nOX.*CpOX +nMX.*CpMX + nB.*CpB + nE.*CpE + nA.*CpA + nW.*CpW;
    dTdv = -QR./Cpchange;
    dPdv = 0;

    %vector of differentials for ODE solver
    dxdv=[dnTdv;dnMdv;dnPXdv;dnOXdv;dnMXdv;dnBdv;dnEdv;dnAdv;dnWdv;dEx1dv;dEx2dv;dEx3dv;dEx4dv;dTdv;dPdv];
end


function dxdv = IsothermalReactor(v,x) 

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

    % %Heat Capacities Data (entered at 6 atm): to be entered: 
    CpT =230.5145;
    CpM =78.6352;
    CpPX =274.1495;
    CpOX =275.6325;
    CpMX=274.1535;
    CpB =187.458;
    CpE =83.01595;
    CpA =88.21465;
    CpW = 38.61675;

    % % Heat of Reaction Data:
    H1 = -69825.5;
    H2 =325.5;
    H3 = 58115.425;
    H4 = 1017.5;
        
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

    %T change
    dTdv = 0;
    dPdv = 0;

    %vector of differentials for ODE solver
    dxdv=[dnTdv;dnMdv;dnPXdv;dnOXdv;dnMXdv;dnBdv;dnEdv;dnAdv;dnWdv;dEx1dv;dEx2dv;dEx3dv;dEx4dv;dTdv;dPdv];
end