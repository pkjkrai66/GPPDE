% Richards Equation
% VanGenuchten
% 10-Feb-2014, 2-June-2016,

function [head, z, deltaZ, deltaT, Error, h] = ...
    Richards_Celia_MixForm_VanGenuchten(Conditions, z, t, thta, Param, Ks)
%{
% Richards_Hform_CeliaFormation solves Richard equation as given by Celia, 1990.
%
% IMPORTANT: Boundary Conditions are defined inside the function.
%-INPUTS ----------------------------
%  Conditions     = Initial and Boundary Conditions as formated by
%                   importConditionsfromXLS.m
%  z.ProfileDepth = Depth of Soil Profile; POSITIVE UPWARD
%  z.discreteZ    = Value of z at which analysis is required
%  z.deltaZ       = grid size for Z, if deltaZ is not a factor of Z then ...
%                   deltaZ is NOT corrected 
%  t.discreteT  = time instanes at which the head is required
%  t.deltaT     = grid size for time, if deltaT is not a factor of T then deltaT 
%                 is NOT corrected
%  thta.r = resudual Soil Moisture
%  thta.s = Saturated Soil Moisture
%  Param   = [alpha N];
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

% [head, z, deltaZ, deltaT, Error, h] = Richards_Celia_MixForm_VanGenuchten(60, 1, 30, .1, .05, 2)
% 
% In followin code variable
%
%}
alpha = Param.alpha;
N = Param.N;
% Ks = 1.5; 
M = 1-1/N;


% Z discritization
Z            = z.discreteZ;
profileDepth = z.ProfileDepth; % Max depth of Soil Profile
deltaZ       = z.deltaZ;

if profileDepth < max(Z)
    error('Values of Vector Z should be less than profileDepth')
end

temp1 = [Z(:); profileDepth(:)];
rmaindrZ = rem(temp1, deltaZ);
rmaindrSumZ = sum(rmaindrZ);


if rmaindrSumZ > 0
    error('deltaZ must be a factor of Vector z')
end

totalzNodes = 1+ profileDepth/deltaZ;
depths = deltaZ*(0:totalzNodes-1);
Z_indices = 1 + round(Z/deltaZ);

% T discritization
T = t.discreteT;
deltaT = t.deltaT;

rmaindrT = rem(T, deltaT);
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

NoIterations = 10000;
h = nan(totalzNodes, totaltimeNodes);     % (i, n) = Distance, Time

topBoundry = interp1q(Conditions.TimeVector, Conditions.TopBoundry, times');
botBoundry = interp1q(Conditions.TimeVector, Conditions.BottomBoundry, times');
initCond = interp1q(Conditions.HeightVect, Conditions.InitiHead, depths');

% NOTE- Height is positive upward.

h(:,1) = initCond;         % Initial Condition
h(1,:) = botBoundry;         % boundtyCdn1 
h(totalzNodes,:) = topBoundry;    % boundtyCdn2

if exist('Conditions.sink', 'var')
    sink = zeros(totalzNodes, totaltimeNodes);
else
    sink = zeros(totalzNodes, totaltimeNodes);
end

ErrorThreshold = .001;
Error = nan(totaltimeNodes,NoIterations);

% Soil Moisture
theta = nan(totalzNodes, totaltimeNodes);
theta(:,1) = SoilMoistureVanGenuchten(h(:,1)); % Fixing initial condition in terms of SM

for n = 1:totaltimeNodes-1
    h(2:totalzNodes-1,n+1) = h(2:totalzNodes-1,n); % Guess value of h@time:n+1 = h@n
    theta(2:totalzNodes-1,n+1) = SoilMoistureVanGenuchten(h(2:totalzNodes-1,n));
    for m = 1:NoIterations
       
        CoeffMat = zeros(totalzNodes, totalzNodes);
        RHSvect = zeros(totalzNodes,1);
        
        for i = 2:(totalzNodes-1)
            
            H1 = h(i-1, n+1);
            H2 = h(i, n+1);
            H3 = h(i+1, n+1);
                   
            K1 = PermeabilityVanGenuchten(H1);
            K2 = PermeabilityVanGenuchten(H2);
            K3 = PermeabilityVanGenuchten(H3);
            
            
            if H2 < 0
				C = -M*(thta.s - thta.r)*(alpha*N*(alpha*abs(H2))^(N-1)*sign(H2))/(1+(alpha*abs(H2))^N)^(M+1);
            else
                C = 0;
            end
            
            CoeffMat(i,i) = C/deltaT + (K1 + 2*K2 + K3)/(2*deltaZ^2);
            RHSvect(i) = C*h(i,n+1)/deltaT + (K3-K1)/(2*deltaZ) - (theta(i,n+1)-theta(i,n))/deltaT - sink(i,n);
            
            
            if (i == 2)
                CoeffMat(i,i+1) = -(K2+K3)/(2*deltaZ^2);
                RHSvect(i) = RHSvect(i) + h(i-1,n+1)*(K1+K2)/(2*deltaZ^2);
            elseif (i == totalzNodes-1)
                CoeffMat(i,i-1) = -(K1+K2)/(2*deltaZ^2);
                RHSvect(i) = RHSvect(i)+ h(i+1,n+1)*(K2+K3)/(2*deltaZ^2);
            else
                CoeffMat(i,i+1) = -(K2+K3)/(2*deltaZ^2);
                CoeffMat(i,i-1) = -(K1+K2)/(2*deltaZ^2);
            end

        end

        X = (CoeffMat(2:totalzNodes-1, 2:totalzNodes-1))\RHSvect(2:totalzNodes-1);
        % ----------If Ceeff. matric is not Positive Definite then
        % ----------use Choleski Decomposition ------------------
%         IsPositiveDefinite = isposdef(CoeffMat(2:zNodes-1, 2:zNodes-1));
%         if IsPositiveDefinite == 1
%             X = (CoeffMat(2:zNodes-1, 2:zNodes-1))\RHSvect(2:zNodes-1);
%         else
%             % --------Using LightSpeed Toolbox ----------     
%             X1 = inv_posdef(CoeffMat(2:zNodes-1, 2:zNodes-1))*RHSvect(2:zNodes-1); %#ok<NASGU>
%             % ------------------
%             fprintf('Matrix is not Positive Definite')
%         end
        
            Error(n,m) = sqrt(sum((h(2:totalzNodes-1, n+1)-X).^2));  % Error(Time, Iteration)
            h(2:totalzNodes-1, n+1) = X;
            theta(1:totalzNodes,n+1) = SoilMoistureVanGenuchten(h(:, n+1));

            if Error(n,m)<ErrorThreshold
                break
            end
    end
    % If code takes many iterations to converge Then Exit the function
    if m == NoIterations
        head = -9999.9999;
        return
    end
    
        
%     fprintf('n = %d   m = %d \n', n, m)
end
head = h(Z_indices, time_Indices);

% -----------------FUNCTIONS----------------

function theta = SoilMoistureVanGenuchten(h)
    tempI = h < 0;
    theta(tempI) = thta.r+ (thta.s-thta.r)*(1./(1+(alpha*abs(h(tempI))).^N)).^M;
    theta(~tempI) = thta.s;
end

function k = PermeabilityVanGenuchten(h)
    tempI = h < 0;
    k(tempI) = Ks*(1-((alpha*abs(h(tempI))).^(N-1)).*((1+(alpha*abs(h(tempI))).^N).^-M)).^2./((1+(alpha*abs(h(tempI))).^N).^(M/2));
    k(~tempI) = Ks;
end

% Movie 

% figure
% axis([-75 -15 0 85])
% box on
% ylabel('Height')
% xlabel('head')


% [rowIndex, columnIndex] = size(h);
% z = 0:1:rowIndex-1;
% h1 = plot(h(:,1), z);
% set(gca, 'xlim', [-75 -15] )
% 
% for i =1:columnIndex
%     set(h1, 'XData', h(:,i), 'YData', z)



%     drawnow
%     pause(.01)
% end

end

