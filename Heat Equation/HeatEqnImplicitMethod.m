% Heat Equation
% 9-Aug-2014

function [head, deltaZ, deltaT, h] = HeatEqnImplicitMethod(Conditions, z, t, k)

% HeatEqnImplicitMethod solves heat equation [Chapra and Canale p. 729].
%
% IMPORTANT: Boundary Conditions are defined inside the function.
%-INPUTS ----------------------------
%  Conditions     = Initial and Boundary Conditions as formated by
%                   importConditionsfromXLS.m
%  z.ProfileDepth = Depth of Soil Profile; POSITIVE UPWARD
%  z.discreteZ    = Value of z at which analysis is required
%  z.deltaZ       = grid size for Z, if deltaZ is not a factor of Z then ...
%                   deltaZ is NOT corrected 
%  t       = time instanes at which the head is required
%  deltaT  = grid size for time, if deltaT is not a factor of T then deltaT 
%            is NOT corrected 
%  k   =    ;
%-OUTPUTS ----------------------------
% head     = head in soil profile at time vector 't' at discrete value of Z
%            Format: h(SPACE, TIME), h(1,:) corrosponds to first coloum gives initial condition.
%                   h(1,:) = Bottom Boundty Cdn
%                   h(zNodes,:) =  Top BoundtyCdn2
% z        = discete values of Z
% deltaZ   = if (input deltaZ is a factor of Z)
%            then deltaZ = deltaZ   
%            else (deltaZ = closest factor of Z)
% h        = Matrix which provides variation of head at all discratized time & space grid points.
%            Format: h(SPACE, TIME), h(:, 1) corrosponds to first coloum gives initial condition.

% [head, z, deltaZ, deltaT, h] = HeatEqnImplicitMethod(60, 1, 30, .1, .05, 2)
% 
% In followin code variable
%


% Z discritization
Z            = z.discreteZ;
barLength = z.length; % Max depth of Soil Profile
deltaZ       = z.deltaZ;
if max(Conditions.HeightVect) > z.length
    error('max(Conditions.HeightVect) can not be more then z.length')
end

if barLength < max(Z)
    error('Values of Vector Z should be less than profileDepth')
end
T = t.discreteT;
deltaT = t.deltaT;
if max(T) > max(Conditions.TimeVector)
    error('maximum value of t must be less then Conditions.TimeVector')
end

temp1 = [Z(:); barLength(:)];
temp11 = round(temp1*10000)/10000;
rmaindrZ = rem(temp11, deltaZ);
rmaindrSumZ = sum(rmaindrZ);

if rmaindrSumZ > 0
    error('deltaZ must be a factor of Vector z')
end

totalzNodes = 1+ round(barLength/deltaZ);
depths = deltaZ*(0:totalzNodes-1);
Z_indices = round((1 + Z/deltaZ)*1000)/1000;

% T discritization
temp2 = round(T*1e4)/1e4;
rmaindrT = rem(temp2, deltaT);
rmaindrSumT = sum(abs(rmaindrT)>eps);

if rmaindrSumT > 0
    error('deltaT must be a factor of Vector t')
end

totaltimeNodes = 1 + round(max(T)/deltaT);
times = deltaT*(0:totaltimeNodes-1);
time_Indices = 1 + round(T/deltaT);

% timeNodes = Time.timeNodes;
% t = Time.t;
% deltaT = T/(timeNodes-1);

h = nan(totalzNodes, totaltimeNodes);     % (i, n) = Distance, Time

topBoundry = interp1q(Conditions.TimeVector', Conditions.TopBoundry', times');
botBoundry = interp1q(Conditions.TimeVector', Conditions.BottomBoundry', times');
initCond = interp1q(Conditions.HeightVect', Conditions.InitiHead', depths');

% NOTE- Height is positive upward.

h(:,1) = initCond;         % Initial Condition
h(1,:) = botBoundry;         % boundtyCdn1 
h(totalzNodes,:) = topBoundry;    % boundtyCdn2

lambda = k*deltaT/(deltaZ^2);

for n = 1:totaltimeNodes-1
    CoeffMat = zeros(totalzNodes, totalzNodes);
    RHSvect = zeros(totalzNodes,1);

    for i = 2:(totalzNodes-1)
        CoeffMat(i, i) = 1+2*lambda;
        RHSvect(i) = h(i,n);
        
        if (i == 2)
            CoeffMat(i,i+1) = -lambda;
            RHSvect(i) = RHSvect(i) + lambda*h(i-1, n+1);
        elseif (i == totalzNodes-1)
            CoeffMat(i,i-1) = -lambda;
            RHSvect(i) = RHSvect(i) + lambda*h(i+1, n+1);
        else
            CoeffMat(i,i+1) = -lambda;
            CoeffMat(i,i-1) = -lambda;
        end

    end
    
    temperature = CoeffMat(2:(totalzNodes-1), 2:(totalzNodes-1))\RHSvect(2:(totalzNodes-1));
    h(2:(totalzNodes-1), n+1) = temperature;
    
end

head = h(Z_indices, time_Indices);

