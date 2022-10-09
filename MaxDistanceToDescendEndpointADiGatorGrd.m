% This code was generated using ADiGator version 1.4
% ©2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function output = MaxDistanceToDescendEndpointADiGatorGrd(input)
global ADiGator_MaxDistanceToDescendEndpointADiGatorGrd
if isempty(ADiGator_MaxDistanceToDescendEndpointADiGatorGrd); ADiGator_LoadData(); end
Gator1Data = ADiGator_MaxDistanceToDescendEndpointADiGatorGrd.MaxDistanceToDescendEndpointADiGatorGrd.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: % tf = input.phase(1).finaltime;
%User Line: % output.objective = tf;
xf.dv = input.phase.finalstate.dv(4);
xf.f = input.phase.finalstate.f(4);
%User Line: xf = input.phase(1).finalstate(4);
output.objective.dv = -xf.dv;
output.objective.f = uminus(xf.f);
%User Line: output.objective = -xf;
%User Line: %----------------------------------%
%User Line: % END: function minimumClimbCost.m %
%User Line: %----------------------------------%
output.objective.dv_size = 10;
output.objective.dv_location = Gator1Data.Index1;
end


function ADiGator_LoadData()
global ADiGator_MaxDistanceToDescendEndpointADiGatorGrd
ADiGator_MaxDistanceToDescendEndpointADiGatorGrd = load('MaxDistanceToDescendEndpointADiGatorGrd.mat');
return
end