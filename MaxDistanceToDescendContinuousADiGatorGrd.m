% This code was generated using ADiGator version 1.4
% �2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function phaseout = MaxDistanceToDescendContinuousADiGatorGrd(input)
global ADiGator_MaxDistanceToDescendContinuousADiGatorGrd
if isempty(ADiGator_MaxDistanceToDescendContinuousADiGatorGrd); ADiGator_LoadData(); end
Gator1Data = ADiGator_MaxDistanceToDescendContinuousADiGatorGrd.MaxDistanceToDescendContinuousADiGatorGrd.Gator1Data;
% ADiGator Start Derivative Computations
g = input.auxdata.g;
%User Line: g    = input.auxdata.g;
m = input.auxdata.m;
%User Line: m    = input.auxdata.m;
S = input.auxdata.S;
%User Line: S    = input.auxdata.S;
rho0 = input.auxdata.rho0;
%User Line: rho0 = input.auxdata.rho0;
a0 = input.auxdata.a0;
%User Line: a0   = input.auxdata.a0;
umax = input.auxdata.umax;
%User Line: umax = input.auxdata.umax;
t.dV = input.phase.time.dV; t.f = input.phase.time.f;
%User Line: t   = input.phase.time;
h.dV = input.phase.state.dV(:,1);
h.f = input.phase.state.f(:,1);
%User Line: h   = input.phase.state(:,1);
v.dV = input.phase.state.dV(:,2);
v.f = input.phase.state.f(:,2);
%User Line: v   = input.phase.state(:,2);
fpa.dV = input.phase.state.dV(:,3);
fpa.f = input.phase.state.f(:,3);
%User Line: fpa = input.phase.state(:,3);
x.dV = input.phase.state.dV(:,4);
x.f = input.phase.state.f(:,4);
%User Line: x   = input.phase.state(:,4);
u.dV = input.phase.control.dV(:,1);
u.f = input.phase.control.f(:,1);
%User Line: u   = input.phase.control(:,1);
throttle.dV = input.phase.control.dV(:,2);
throttle.f = input.phase.control.f(:,2);
%User Line: throttle = input.phase.control(:,2);
hbar.dV = h.dV./1000;
hbar.f = h.f/1000;
%User Line: hbar = h./1000;
cada1f1 = size(t.f,1);
HH.f = zeros(cada1f1,5);
%User Line: HH   = zeros(size(t,1),5);
cada1f1 = size(t.f);
cada1f2 = ones(cada1f1);
HH.f(:,1) = cada1f2;
%User Line: HH(:,1) = ones(size(t));
cadaforvar1.f =  2:6;
%User Line: cadaforvar1 = 2:6;
HH.dV = zeros(size(HH.f,1),5);
HH.f(:,6) = 0;
for cadaforcount1 = 1:5
    ii.f = cadaforvar1.f(:,cadaforcount1);
    %User Line: ii = cadaforvar1(:,cadaforcount1);
    cada1f1 = ii.f - 1;
    cada1td1 = zeros(size(HH.f,1),1);
    cada1td1(:,logical(Gator1Data.Index1(:,cadaforcount1))) = HH.dV(:,nonzeros(Gator1Data.Index1(:,cadaforcount1)));
    cada1f2dV = cada1td1;
    cada1f2 = HH.f(:,cada1f1);
    cada1td1 = hbar.f.*cada1f2dV;
    cada1td1 = cada1td1 + cada1f2.*hbar.dV;
    cada1f3dV = cada1td1;
    cada1f3 = cada1f2.*hbar.f;
    HH.dV(:,logical(Gator1Data.Index2(:,cadaforcount1))) = cada1f3dV(:,nonzeros(Gator1Data.Index2(:,cadaforcount1)));
    HH.f(:,ii.f) = cada1f3;
    %User Line: HH(:,ii) = HH(:,ii-1).*hbar;
end
cada1f1dV = HH.dV(:,Gator1Data.Index6);
cada1f1 = HH.f(:,Gator1Data.Index5);
cada1tf1 = zeros(4,1);
cada1tf1(Gator1Data.Index8) = input.auxdata.z(Gator1Data.Index7);
z.dV = cada1f1dV*cada1tf1;
z.f = cada1f1*input.auxdata.z;
%User Line: z = HH(:,2:5)*input.auxdata.z;
cada1f1dV = -z.dV;
cada1f1 = uminus(z.f);
cada1f2dV = exp(cada1f1).*cada1f1dV;
cada1f2 = exp(cada1f1);
r.dV = input.auxdata.r0.*cada1f2dV;
r.f = input.auxdata.r0*cada1f2;
%User Line: r = input.auxdata.r0*exp(-z);
cada1f1dV = HH.dV(:,1);
cada1f1 = HH.f(:,Gator1Data.Index9);
cada1tf1 = zeros(1,1);
cada1tf1(1) = input.auxdata.y(2);
cada1f2dV = cada1f1dV*cada1tf1;
cada1f2 = cada1f1*input.auxdata.y;
cada1td1 = cada1f2dV;
cada1td1 = cada1td1 + r.dV;
y.dV = cada1td1;
y.f = cada1f2 + r.f;
%User Line: y = HH(:,1:2)*input.auxdata.y + r;
cada1f1dV = exp(y.f).*y.dV;
cada1f1 = exp(y.f);
rho.dV = rho0.*cada1f1dV;
rho.f = rho0*cada1f1;
%User Line: rho = rho0*exp(y);
cada1f1dV = HH.dV(:,Gator1Data.Index11);
cada1f1 = HH.f(:,Gator1Data.Index10);
cada1tf1 = zeros(3,1);
cada1tf1(Gator1Data.Index13) = input.auxdata.theta(Gator1Data.Index12);
theta.dV = cada1f1dV*cada1tf1;
theta.f = cada1f1*input.auxdata.theta;
%User Line: theta = HH(:,1:4)*input.auxdata.theta;
cada1f1dV = (1/2)./sqrt(theta.f).*theta.dV;
cada1f1dV(theta.f == 0 & theta.dV == 0) = 0;
cada1f1 = sqrt(theta.f);
a.dV = a0.*cada1f1dV;
a.f = a0*cada1f1;
%User Line: a     = a0*sqrt(theta);
cada1td1 = zeros(size(v.dV,1),2);
cada1td1(:,2) = v.dV./a.f;
cada1td1(:,1) = cada1td1(:,1) + -v.f./a.f.^2.*a.dV;
M.dV = cada1td1;
M.f = v.f./a.f;
%User Line: M = v./a;
cada1td1 = zeros(size(rho.dV,1),2);
cada1td1(:,1) = v.f.*rho.dV;
cada1td1(:,2) = cada1td1(:,2) + rho.f.*v.dV;
cada1f1dV = cada1td1;
cada1f1 = rho.f.*v.f;
cada1tf1 = v.f(:,Gator1Data.Index14);
cada1td1 = cada1tf1.*cada1f1dV;
cada1td1(:,2) = cada1td1(:,2) + cada1f1.*v.dV;
cada1f2dV = cada1td1;
cada1f2 = cada1f1.*v.f;
cada1f3dV = input.auxdata.S.*cada1f2dV;
cada1f3 = cada1f2*input.auxdata.S;
q.dV = cada1f3dV./2;
q.f = cada1f3/2;
%User Line: q = rho.*v.*v.*input.auxdata.S/2;
cada1f1 = m*g;
L.dV = cada1f1.*u.dV;
L.f = cada1f1*u.f;
%User Line: L = m*g*u;
cada1f1 = size(t.f,1);
MM.f = zeros(cada1f1,6);
%User Line: MM   = zeros(size(t,1),6);
cada1f1 = size(t.f);
cada1f2 = ones(cada1f1);
MM.f(:,1) = cada1f2;
%User Line: MM(:,1) = ones(size(t));
cadaforvar2.f =  2:6;
%User Line: cadaforvar2 = 2:6;
MM.dV = zeros(size(MM.f,1),10);
for cadaforcount2 = 1:5
    ii.f = cadaforvar2.f(:,cadaforcount2);
    %User Line: ii = cadaforvar2(:,cadaforcount2);
    cada1f1 = ii.f - 1;
    cada1td1 = zeros(size(MM.f,1),2);
    cada1td1(:,logical(Gator1Data.Index3(:,cadaforcount2))) = MM.dV(:,nonzeros(Gator1Data.Index3(:,cadaforcount2)));
    cada1f2dV = cada1td1;
    cada1f2 = MM.f(:,cada1f1);
    cada1tf1 = M.f(:,Gator1Data.Index15);
    cada1td1 = cada1tf1.*cada1f2dV;
    cada1tf1 = cada1f2(:,Gator1Data.Index16);
    cada1td1 = cada1td1 + cada1tf1.*M.dV;
    cada1f3dV = cada1td1;
    cada1f3 = cada1f2.*M.f;
    MM.dV(:,logical(Gator1Data.Index4(:,cadaforcount2))) = cada1f3dV(:,nonzeros(Gator1Data.Index4(:,cadaforcount2)));
    MM.f(:,ii.f) = cada1f3;
    %User Line: MM(:,ii) = MM(:,ii-1).*M;
end
cada1f1dV = MM.dV(:,Gator1Data.Index18);
cada1f1 = MM.f(:,Gator1Data.Index17);
cada1tf1 = zeros(10,2);
cada1tf1(Gator1Data.Index20) = input.auxdata.a(Gator1Data.Index19);
numeratorCD0.dV = cada1f1dV*cada1tf1;
numeratorCD0.f = cada1f1*input.auxdata.a;
%User Line: numeratorCD0   = MM(:,1:6)*input.auxdata.a;
cada1f1dV = MM.dV(:,Gator1Data.Index22);
cada1f1 = MM.f(:,Gator1Data.Index21);
cada1tf1 = zeros(10,2);
cada1tf1(Gator1Data.Index24) = input.auxdata.b(Gator1Data.Index23);
denominatorCD0.dV = cada1f1dV*cada1tf1;
denominatorCD0.f = cada1f1*input.auxdata.b;
%User Line: denominatorCD0 = MM(:,1:6)*input.auxdata.b;
cada1tf1 = denominatorCD0.f(:,Gator1Data.Index25);
cada1td1 = numeratorCD0.dV./cada1tf1;
cada1tf1 = numeratorCD0.f(:,Gator1Data.Index26);
cada1tf2 = denominatorCD0.f(:,Gator1Data.Index27);
cada1td1 = cada1td1 + -cada1tf1./cada1tf2.^2.*denominatorCD0.dV;
CD0.dV = cada1td1;
CD0.f = numeratorCD0.f./denominatorCD0.f;
%User Line: CD0            = numeratorCD0./denominatorCD0;
cada1f1dV = MM.dV(:,Gator1Data.Index29);
cada1f1 = MM.f(:,Gator1Data.Index28);
cada1tf1 = zeros(10,2);
cada1tf1(Gator1Data.Index31) = input.auxdata.c(Gator1Data.Index30);
numeratorK.dV = cada1f1dV*cada1tf1;
numeratorK.f = cada1f1*input.auxdata.c;
%User Line: numeratorK     = MM(:,1:6)*input.auxdata.c;
cada1f1dV = MM.dV(:,Gator1Data.Index33);
cada1f1 = MM.f(:,Gator1Data.Index32);
cada1tf1 = zeros(10,2);
cada1tf1(Gator1Data.Index35) = input.auxdata.d(Gator1Data.Index34);
denominatorK.dV = cada1f1dV*cada1tf1;
denominatorK.f = cada1f1*input.auxdata.d;
%User Line: denominatorK   = MM(:,1:6)*input.auxdata.d;
cada1tf1 = denominatorK.f(:,Gator1Data.Index36);
cada1td1 = numeratorK.dV./cada1tf1;
cada1tf1 = numeratorK.f(:,Gator1Data.Index37);
cada1tf2 = denominatorK.f(:,Gator1Data.Index38);
cada1td1 = cada1td1 + -cada1tf1./cada1tf2.^2.*denominatorK.dV;
K.dV = cada1td1;
K.f = numeratorK.f./denominatorK.f;
%User Line: K              = numeratorK./denominatorK;
cada1f1 = m^2;
cada1f2 = g^2;
cada1f3 = cada1f1*cada1f2;
cada1tf2 = q.f(:,Gator1Data.Index39);
cada1f4dV = 2.*cada1tf2.^(2-1).*q.dV;
cada1f4 = q.f.^2;
cada1tf2 = cada1f4(:,Gator1Data.Index40);
cada1f5dV = -cada1f3./cada1tf2.^2.*cada1f4dV;
cada1f5 = cada1f3./cada1f4;
cada1tf1 = cada1f5(:,Gator1Data.Index41);
cada1td1 = cada1tf1.*K.dV;
cada1tf1 = K.f(:,Gator1Data.Index42);
cada1td1 = cada1td1 + cada1tf1.*cada1f5dV;
cada1f6dV = cada1td1;
cada1f6 = K.f.*cada1f5;
cada1f7dV = 2.*u.f.^(2-1).*u.dV;
cada1f7 = u.f.^2;
cada1tf1 = cada1f7(:,Gator1Data.Index43);
cada1td1 = zeros(size(cada1f6dV,1),3);
cada1td1(:,Gator1Data.Index44) = cada1tf1.*cada1f6dV;
cada1td1(:,3) = cada1td1(:,3) + cada1f6.*cada1f7dV;
cada1f8dV = cada1td1;
cada1f8 = cada1f6.*cada1f7;
cada1td1 = zeros(size(CD0.dV,1),3);
cada1td1(:,Gator1Data.Index45) = CD0.dV;
cada1td1 = cada1td1 + cada1f8dV;
cada1f9dV = cada1td1;
cada1f9 = CD0.f + cada1f8;
cada1tf1 = cada1f9(:,Gator1Data.Index46);
cada1td1 = zeros(size(q.dV,1),3);
cada1td1(:,Gator1Data.Index47) = cada1tf1.*q.dV;
cada1tf1 = q.f(:,Gator1Data.Index48);
cada1td1 = cada1td1 + cada1tf1.*cada1f9dV;
D.dV = cada1td1;
D.f = q.f.*cada1f9;
%User Line: D = q.*(CD0+K.*((m.^2).*(g.^2)./(q.^2)).*(u.^2));
cada1f1 = m^2;
cada1f2 = g^2;
cada1f3 = cada1f1*cada1f2;
cada1tf2 = q.f(:,Gator1Data.Index49);
cada1f4dV = 2.*cada1tf2.^(2-1).*q.dV;
cada1f4 = q.f.^2;
cada1tf2 = cada1f4(:,Gator1Data.Index50);
cada1f5dV = -cada1f3./cada1tf2.^2.*cada1f4dV;
cada1f5 = cada1f3./cada1f4;
cada1tf1 = cada1f5(:,Gator1Data.Index51);
cada1td1 = cada1tf1.*K.dV;
cada1tf1 = K.f(:,Gator1Data.Index52);
cada1td1 = cada1td1 + cada1tf1.*cada1f5dV;
cada1f6dV = cada1td1;
cada1f6 = K.f.*cada1f5;
cada1f7 = umax^2;
cada1f8dV = cada1f7.*cada1f6dV;
cada1f8 = cada1f6*cada1f7;
cada1td1 = CD0.dV;
cada1td1 = cada1td1 + cada1f8dV;
cada1f9dV = cada1td1;
cada1f9 = CD0.f + cada1f8;
cada1tf1 = cada1f9(:,Gator1Data.Index53);
cada1td1 = cada1tf1.*q.dV;
cada1tf1 = q.f(:,Gator1Data.Index54);
cada1td1 = cada1td1 + cada1tf1.*cada1f9dV;
Dmax.dV = cada1td1;
Dmax.f = q.f.*cada1f9;
%User Line: Dmax = q.*(CD0+K.*((m.^2).*(g.^2)./(q.^2)).*(umax.^2));
cada1f1dV = MM.dV(:,Gator1Data.Index56);
cada1f1 = MM.f(:,Gator1Data.Index55);
cada1tf1 = zeros(10,12);
cada1tf1(Gator1Data.Index58) = input.auxdata.f(Gator1Data.Index57);
E.dV = cada1f1dV*cada1tf1;
E.f = cada1f1*input.auxdata.f;
%User Line: E = MM(:,1:6)*input.auxdata.f;
%User Line: % T  = sum(HH(:,1:6).*E,2).*g/2.2;  (original comment)
cada1f1dV = HH.dV(:,Gator1Data.Index60);
cada1f1 = HH.f(:,Gator1Data.Index59);
cada1tf1 = E.f(:,Gator1Data.Index61);
cada1td1 = zeros(size(cada1f1dV,1),12);
cada1td1(:,Gator1Data.Index62) = cada1tf1.*cada1f1dV;
cada1tf1 = cada1f1(:,Gator1Data.Index63);
cada1td1 = cada1td1 + cada1tf1.*E.dV;
cada1f2dV = cada1td1;
cada1f2 = cada1f1.*E.f;
cada1tf2 = ones(6,1);
cada1tf1 = zeros(12,2);
cada1tf1(Gator1Data.Index65) = cada1tf2(Gator1Data.Index64);
T.dV = cada1f2dV*cada1tf1;
T.f = sum(cada1f2,2);
%User Line: T  = sum(HH(:,1:6).*E,2);
cada1f1dV = cos(fpa.f).*fpa.dV;
cada1f1 = sin(fpa.f);
cada1td1 = zeros(size(v.dV,1),2);
cada1td1(:,1) = cada1f1.*v.dV;
cada1td1(:,2) = cada1td1(:,2) + v.f.*cada1f1dV;
hdot.dV = cada1td1;
hdot.f = v.f.*cada1f1;
%User Line: hdot = v.*sin(fpa);
%User Line: %vdot = (throttle.*T-D)./m-g.*sin(fpa);
cada1td1 = zeros(size(T.dV,1),3);
cada1td1(:,Gator1Data.Index66) = T.dV;
cada1td1 = cada1td1 + -D.dV;
cada1f1dV = cada1td1;
cada1f1 = T.f - D.f;
cada1td1 = cada1f1dV;
cada1td1(:,Gator1Data.Index67) = cada1td1(:,Gator1Data.Index67) + Dmax.dV;
cada1f2dV = cada1td1;
cada1f2 = cada1f1 + Dmax.f;
cada1td1 = zeros(size(throttle.dV,1),4);
cada1td1(:,4) = cada1f2.*throttle.dV;
cada1tf1 = throttle.f(:,Gator1Data.Index68);
cada1td1(:,Gator1Data.Index69) = cada1td1(:,Gator1Data.Index69) + cada1tf1.*cada1f2dV;
cada1f3dV = cada1td1;
cada1f3 = throttle.f.*cada1f2;
cada1td1 = cada1f3dV;
cada1td1(:,Gator1Data.Index70) = cada1td1(:,Gator1Data.Index70) + -Dmax.dV;
cada1f4dV = cada1td1;
cada1f4 = cada1f3 - Dmax.f;
cada1f5dV = cada1f4dV./m;
cada1f5 = cada1f4/m;
cada1f6dV = cos(fpa.f).*fpa.dV;
cada1f6 = sin(fpa.f);
cada1f7dV = g.*cada1f6dV;
cada1f7 = g*cada1f6;
cada1td1 = zeros(size(cada1f5dV,1),5);
cada1td1(:,Gator1Data.Index71) = cada1f5dV;
cada1td1(:,3) = cada1td1(:,3) + -cada1f7dV;
vdot.dV = cada1td1;
vdot.f = cada1f5 - cada1f7;
%User Line: vdot = (throttle.*(T-D+Dmax)-Dmax)./m-g.*sin(fpa);
cada1f1dV = -sin(fpa.f).*fpa.dV;
cada1f1 = cos(fpa.f);
cada1td1 = zeros(size(u.dV,1),2);
cada1td1(:,2) = u.dV;
cada1td1(:,1) = cada1td1(:,1) + -cada1f1dV;
cada1f2dV = cada1td1;
cada1f2 = u.f - cada1f1;
cada1f3dV = g.*cada1f2dV;
cada1f3 = g*cada1f2;
cada1tf1 = v.f(:,Gator1Data.Index72);
cada1td1 = zeros(size(cada1f3dV,1),3);
cada1td1(:,Gator1Data.Index73) = cada1f3dV./cada1tf1;
cada1td1(:,1) = cada1td1(:,1) + -cada1f3./v.f.^2.*v.dV;
fpadot.dV = cada1td1;
fpadot.f = cada1f3./v.f;
%User Line: fpadot = g.*(u-cos(fpa))./v;
cada1f1dV = -sin(fpa.f).*fpa.dV;
cada1f1 = cos(fpa.f);
cada1td1 = zeros(size(v.dV,1),2);
cada1td1(:,1) = cada1f1.*v.dV;
cada1td1(:,2) = cada1td1(:,2) + v.f.*cada1f1dV;
xdot.dV = cada1td1;
xdot.f = v.f.*cada1f1;
%User Line: xdot = v.*cos(fpa);
%User Line: % input.auxdata.realthrottle = (vdot.*(m-g.*sin(fpa))+D)./T;
%User Line: %phaseout(1).dynamics = [hdot, vdot, fpadot];
cada1td1 = zeros(size(hdot.f,1),12);
cada1td1(:,Gator1Data.Index74) = hdot.dV;
cada1td1(:,Gator1Data.Index75) = vdot.dV;
cada1td1(:,Gator1Data.Index76) = fpadot.dV;
cada1td1(:,Gator1Data.Index77) = xdot.dV;
cada1f1dV = cada1td1;
cada1f1 = [hdot.f vdot.f fpadot.f xdot.f];
phaseout(1).dynamics.dV = cada1f1dV;
phaseout(1).dynamics.f = cada1f1;
%User Line: phaseout(1).dynamics = [hdot, vdot, fpadot, xdot];
phaseout(1).path.dV = q.dV;
phaseout(1).path.f = q.f;
%User Line: phaseout(1).path = q;
cada1td1 = zeros(size(throttle.dV,1),3);
cada1td1(:,3) = T.f.*throttle.dV;
cada1tf1 = throttle.f(:,Gator1Data.Index78);
cada1td1(:,Gator1Data.Index79) = cada1td1(:,Gator1Data.Index79) + cada1tf1.*T.dV;
cada1f1dV = cada1td1;
cada1f1 = throttle.f.*T.f;
cada1f2dV = -throttle.dV;
cada1f2 = 1 - throttle.f;
cada1td1 = D.dV;
cada1td1(:,Gator1Data.Index80) = cada1td1(:,Gator1Data.Index80) + -Dmax.dV;
cada1f3dV = cada1td1;
cada1f3 = D.f - Dmax.f;
cada1td1 = zeros(size(cada1f2dV,1),4);
cada1td1(:,4) = cada1f3.*cada1f2dV;
cada1tf1 = cada1f2(:,Gator1Data.Index81);
cada1td1(:,Gator1Data.Index82) = cada1td1(:,Gator1Data.Index82) + cada1tf1.*cada1f3dV;
cada1f4dV = cada1td1;
cada1f4 = cada1f2.*cada1f3;
cada1td1 = zeros(size(cada1f1dV,1),4);
cada1td1(:,Gator1Data.Index83) = cada1f1dV;
cada1td1 = cada1td1 + cada1f4dV;
cada1f5dV = cada1td1;
cada1f5 = cada1f1 + cada1f4;
cada1tf1 = T.f(:,Gator1Data.Index84);
cada1td1 = cada1f5dV./cada1tf1;
cada1tf1 = cada1f5(:,Gator1Data.Index85);
cada1tf2 = T.f(:,Gator1Data.Index86);
cada1td1(:,Gator1Data.Index87) = cada1td1(:,Gator1Data.Index87) + -cada1tf1./cada1tf2.^2.*T.dV;
cada1f6dV = cada1td1;
cada1f6 = cada1f5./T.f;
phaseout(1).realthrottle.dV = cada1f6dV;
phaseout(1).realthrottle.f = cada1f6;
%User Line: phaseout(1).realthrottle = (throttle.*T + (1-throttle).*(D-Dmax))./T;
%User Line: %---------------------------------%
%User Line: % END: function maxDistanceDae.m %
%User Line: %---------------------------------%
phaseout.dynamics.dV_size = [4,7];
phaseout.dynamics.dV_location = Gator1Data.Index88;
phaseout.path.dV_size = 7;
phaseout.path.dV_location = Gator1Data.Index89;
phaseout.realthrottle.dV_size = 7;
phaseout.realthrottle.dV_location = Gator1Data.Index90;
end


function ADiGator_LoadData()
global ADiGator_MaxDistanceToDescendContinuousADiGatorGrd
ADiGator_MaxDistanceToDescendContinuousADiGatorGrd = load('MaxDistanceToDescendContinuousADiGatorGrd.mat');
return
end