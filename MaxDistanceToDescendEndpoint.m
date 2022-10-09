function output = MaxDistanceToDescendEndpoint(input)

% tf = input.phase(1).finaltime;
% output.objective = tf;
xf = input.phase(1).finalstate(4);
output.objective = -xf;
%----------------------------------%
% END: function minimumClimbCost.m %
%----------------------------------%
