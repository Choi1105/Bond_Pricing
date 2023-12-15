%% Parameter constraint
% input: psi, Sn (paramters, structure parameters)
% output: valid = 1, parameter constraint satisfied
%         valid = 0, otherwise

function [valid] = paramconst(psi,Sn)
% Define structural parameters
indR = Sn.indR;
indF = Sn.indF;
indQ = Sn.indQ;

% Pre-allocation
validm = ones(30,1);

% Except for "NAN" parameters 
if maxc(isnan(psi)') == 1
    validm(1) = 0;
end 

% variance > 0
validm(2) = minc(psi(indR)) > -5;
validm(3) = minc(psi(indQ)) > -5;

% Stationarity condition
validm(4) = abs(psi(indF)) < 1;  

% Check whether all constranit satisfied
valid = minc(validm); 

end