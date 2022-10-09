% -------------------------------------------------------------------%
% Minimum Time-to-Climb Problem Using Seywald Model.                 %
% The vehicle model for this problem is taken from the               % 
% following two references:                                          %
%   Seywald, H., Cliff, E. M., and Well, K. H.,                      %
%   "Range Optimal Trajectories for an Aircraft Flying in            %
%    the Vertical Plane," Journal of Guidance, Control, and Dynamics,%
%    Vol. 17, No. 2, March-April, 1994.                              %
% and                                                                %
%    Rao, A. V., Extension of a Computational Singular Perturbation  %
%    Methodology to Optimal Control Problems, Ph.D. Thesis, Dept.    %
%    of Mechanical and Aerospace Engineering, Princeton University,  %
%    June 1996.                                                      %
% -------------------------------------------------------------------%
clear all; clc

% -------------------------------------------------------------------%
% Below is a listing of the coefficients used in the model developed %
% by Seywald, Cliff, and Well. It is important to note that MATLAB   %
% indexing starts at unity (not zero).  Some of the coefficients     %
% shown below start with an index of zero in the mathematical model  %
% (to match the power of the quantity by which the coefficient is    %
% multiplied).  In the case of a zero index in the model, the MATLAB % 
% code needs to be shifted appropriately so that the "zero index"    %
% in the model corresponds to an index of "one" in the MATLAB code.  %  
% -------------------------------------------------------------------%
% -------------------------------------------------------------------%
% a_j, (j=0,...4)
% -------------------------------------------------------------------%
a      = [2.61059846050e-2; -8.57043966269e-2; 1.07863115049e-1; -6.44772018636e-2; 1.64933626507e-2; 0];
% -------------------------------------------------------------------%
% b_j, (j=0,...4)                                                    %
% -------------------------------------------------------------------%
b      = [1.37368651246e0; -4.57116286752e0; 5.72789877344e0; -3.25219000620e0; 7.29821847445e-1; 0];
% -------------------------------------------------------------------%
% c_j, (j=0,...4)                                                    %
% -------------------------------------------------------------------%
c      = [1.23001735612e0; -2.97244144190e0; 2.78009092756e0; -1.16227834301e0; 1.81868987624e-1; 0];
% -------------------------------------------------------------------%
% d_j, (j=0,...4)                                                    %
% -------------------------------------------------------------------%
d      = [1.42392902737e1; -3.24759126471e1; 2.96838643792e1; -1.33316812491e1; 2.87165882405e0; -2.27239723756e-1];
% -------------------------------------------------------------------%
% f_{ji}, (i,j=0,...,5).  Each row below is the jth row f(j,:)       %
% -------------------------------------------------------------------%
f(1,:) =  [+533570.719821932 -65279.5544959089 -2029.73573491603 +2208.4885382122  -206.176710590188 +5.3492958964598];
f(2,:) =  [-1569835.98929465 +230941.307799122 +10316.5819842132 -10021.6431737966 +931.394759776983 -23.9850228440534];
f(3,:) =  [+2694697.21158 -426130.648897193 -17322.254298181  +17728.6056742698 -1641.98911534889 +42.1370748538242];
f(4,:) =  [-1918670.43511233 +371189.844582602 +5508.27423299061 -13699.9755247614 +1310.0289617555  -33.9687773600601];
f(5,:) =  [+608767.318797674 -146511.01106585  +2477.19221622313 +4740.84430757305 -480.745025774656 +12.7275592676088];
f(6,:) =  [-74209.5599831021 +21887.7903912124 -1051.60186310807 -607.419609387082 +66.3287011203438 -1.79995568442059];
% -------------------------------------------------------------------%
% y_j, (j=0,1)                                                       %
% -------------------------------------------------------------------%
y      = [-1.0228055; -0.12122693];
% -------------------------------------------------------------------%
% z_j, (j=1,...,4)                                                   %
% -------------------------------------------------------------------%
z      = [-3.48643241e-2; 3.50991865e-3; -8.33000535e-5; 1.15219733e-6];
% -------------------------------------------------------------------%
% theta_j, (j=0,...,3)                                               %
% -------------------------------------------------------------------%
theta  = [292.1; -8.87743; 0.193315; 3.72e-3];


%--------------------------------------------------------------------%
%------------- Physical Constants and Other Parameters --------------%
%--------------------------------------------------------------------%
g = 9.80665;  % Acceleration Due to Gravity (m/s)
m = 16818.18; % Vehicle Mass (kg) 
S = 60;       % Vehicle Reference Area (m^2)
a0 = 20.0468; % Speed of Sound Multiplier (m/s)

%--------------------------------------------------------------------%
%---- All Data Required by User Functions is Stored in AUXDATA ------%
%--------------------------------------------------------------------%
auxdata.a     = a;
auxdata.b     = b;
auxdata.c     = c;
auxdata.d     = d;
auxdata.f     = f;
auxdata.y     = y;
auxdata.z     = z;
auxdata.theta = theta;
auxdata.r0    = 1.0228066;
auxdata.rho0  = 1.225;
auxdata.g     = g;
auxdata.m     = m;
auxdata.S     = S;
auxdata.a0    = a0;

%-------------------------------------------------------------------%
%----------------------- Boundary Conditions -----------------------%
%-------------------------------------------------------------------%
t0 = 0;           % Initial Time (sec)
%tf = 1000;       % Final Time (sec)
tf = 60;
%h0 = 0;        % Initial Altitude (m)
h0 = 12119.324;
%hf = 19994.88/2; % Final Altitude (m)
hf = 942.292;
% hf = 10000;
% hf = 12119;
% hf = 500;
%v0 = 129.314;  % Initial Speed (m/s)
v0 = 712.866;
%vf = 295.092;  % Terminal Speed (m/s)
vf = 397.54;
fpa0 = 0;      % Initial Flight Path Angle (radians)
%fpaf = 0;      % Terminal Flight Path Angle (radians)
fpaf = -0.2;
x0 = 0;        % Initial Position (m)

%-------------------------------------------------------------------%
%----------------------- Limits on Variables -----------------------%
%-------------------------------------------------------------------%
hmin   =  0;          hmax   =  21031.2; 
%vmin   =  3.048;      vmax   =  609.6;
vmin   = 10;          vmax   = 2000;
fpamin = -90*pi/180;  fpamax = -fpamin;
xmin = 0;             xmax   = 100000;
umin   = -10;         umax   =  10;
tfmin  = 10;          tfmax  = 1000;
qmin   = 0;           qmax   = 96535.5*S;
throttlemin = 0;      throttlemax = 1;

auxdata.umax = umax;
% Phase 1 Information
iphase = 1;
bounds.phase(iphase).initialtime.lower  = t0;
bounds.phase(iphase).initialtime.upper  = t0;
bounds.phase(iphase).finaltime.lower    = tf;
bounds.phase(iphase).finaltime.upper    = tf;
% bounds.phase(iphase).initialstate.lower = [h0, v0, fpa0];
% bounds.phase(iphase).initialstate.upper = [h0, v0, fpa0];
% bounds.phase(iphase).state.lower        = [hmin, vmin, fpamin];
% bounds.phase(iphase).state.upper        = [hmax, vmax, fpamax];
% bounds.phase(iphase).finalstate.lower   = [hf, vf, fpaf];
% bounds.phase(iphase).finalstate.upper   = [hf, vf, fpaf];
bounds.phase(iphase).initialstate.lower = [h0, v0, fpa0, x0];
bounds.phase(iphase).initialstate.upper = [h0, v0, fpa0, x0];
bounds.phase(iphase).state.lower        = [hmin, vmin, fpamin, xmin];
bounds.phase(iphase).state.upper        = [hmax, vmax, fpamax, xmax];
bounds.phase(iphase).finalstate.lower   = [hf, vf, fpaf, xmin];
bounds.phase(iphase).finalstate.upper   = [hf, vf, fpaf, xmax];
bounds.phase(iphase).path.lower        = [qmin];
bounds.phase(iphase).path.upper        = [qmax];

bounds.phase(iphase).control.lower     = [umin throttlemin];
bounds.phase(iphase).control.upper     = [umax throttlemax];
guess.phase(iphase).time                = [t0; tfmax];
guess.phase(iphase).state(:,1)          = [h0; hf];
guess.phase(iphase).state(:,2)          = [v0; vf];
guess.phase(iphase).state(:,3)          = [fpa0; fpaf];
guess.phase(iphase).state(:,4)          = [x0; xmax];
guess.phase(iphase).control(:,1)        = [-1; 1];
guess.phase(iphase).control(:,2)        = [1; 1];

%-------------------------------------------------------------------------%
%----------Provide Mesh Refinement Method and Initial Mesh ---------------%
%-------------------------------------------------------------------------%
mesh.method       = 'hp-PattersonRao';
mesh.tolerance    = 1e-6;
mesh.maxiterations = 5;
mesh.colpointsmin = 4;
mesh.colpointsmax = 10;
NN = 10;
mesh.phase.colpoints = 4*ones(1,NN);
mesh.phase.fraction = ones(1,NN)/NN;

%-------------------------------------------------------------------%
%---------- Configure Setup Using the information provided ---------%
%-------------------------------------------------------------------%
setup.name                             = 'Max-Distance-To-Descend-Problem';
setup.functions.continuous             = 'MaxDistanceToDescendContinuous';
setup.functions.endpoint               = 'MaxDistanceToDescendEndpoint';
setup.nlp.solver                       = 'ipopt';
setup.nlp.ipoptoptions.tolerance       = 1e-7;
setup.nlp.ipoptoptions.linear_solver   = 'ma57';
setup.displaylevel                     = 2;
setup.bounds                           = bounds;
setup.guess                            = guess;
setup.mesh                             = mesh;
setup.auxdata                          = auxdata;
setup.derivatives.supplier             = 'adigator';
setup.derivatives.derivativelevel      = 'second';
setup.scales.method                    = 'none';
setup.method                           = 'RPM-Differentiation';

output = gpops2(setup);

%-----------------------------------------%
% END: function seywaldMinimumClimbMain.m %
%-----------------------------------------%

%Plot the output states versus time

%Create dynamic pressure graph from paper(0.5*density*v^2)
%Velocity versus height
hdata = output.result.solution.phase.state(:,1);
vdata = output.result.solution.phase.state(:,2);

%Create dynamic pressure limit versus altitude (0.5*density*v^2)
h = h0+500:-10:hmin;

%----------------------------------preliminary constants------------------%

z =  -3.48643241*10^(-5).*h + 3.50991865*10^(-9).*h.^2 + -8.33000535*10^(-14).*h.^3 + 1.15219733*10^(-18).*h.^4;
r =  1.0228055.*exp(-z);
y = -0.12122693.*10^(-3).*h + r - 1.0228055;
%-------------------------------------------------------------------------%

%Compute density at altitude
density = 1.225 .* exp(y);
v = sqrt(2.*qmax ./ (S.*density));

figure
plot(vdata,hdata,'LineWidth',3);
hold on
plot(v,h);
title('altitude versus velocity')
legend('flight path','dynamic pressure limit')
xlabel('V (m/s)')
ylabel('Altitude (m)')

figure
pp = plot(output.result.solution.phase.time,output.result.solution.phase(1).control);
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u(t)$','Interpreter','LaTeX');
set(xl,'Fontsize',18);
set(yl,'Fontsize',18);
set(gca,'Fontsize',16);
set(pp,'LineWidth',1.25);
title('Control Vs. Time')
grid on

figure
plot(output.result.solution.phase.state(:,4),output.result.solution.phase.state(:,1));
title('Altitude versus Downrange')
xlabel('Downrange (m)')
ylabel('Altitude (m)')

figure
ax1 = subplot(221);
ph= plot(output.result.solution.phase.time,output.result.solution.phase.state(:,1));
set(ph,'linewidth',2)
set(gca,'fontsize',18)
xlabel('Time (s)')
ylabel('Altitude (m)');
th = title('States Versus Time');
set(th,'fontsize',20,'fontweight','bold')

subplot(222)
ph= plot(output.result.solution.phase.time,output.result.solution.phase.state(:,2));
set(ph,'linewidth',2)
set(gca,'fontsize',18)
xlabel('Time (s)')
ylabel('Velocity (m/s)');

subplot(223)
ph= plot(output.result.solution.phase.time,output.result.solution.phase.state(:,3));
set(ph,'linewidth',2)
set(gca,'fontsize',18)
xlabel('Time (s)');
ylabel('Flight Path Angle (radians)');


subplot(224)
ph= plot(output.result.solution.phase.time,output.result.solution.phase.state(:,4));
set(ph,'linewidth',2)
set(gca,'fontsize',18)
xlabel('Time (s)');
ylabel('Downrange (m)');
 



