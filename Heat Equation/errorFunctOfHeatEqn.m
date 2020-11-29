function [errNorm2, errNorm1, errNormInf] = errorFunctOfHeatEqn(eqnParam,...
                      delY_delt_GPMat, del2Y_delz2_GPMat,  VarianceMatrix, DataTrim)

% [errNorm2, errNorm1, errNormInf] = errorFunctOfHeatEqn(eqnParam,...
%                       delY_delt_GPMat, del2Y_delz2_GPMat,  VarianceMatrix, DataTrim)

% if (nargin <4)|| isempty(VarianceMatrix)
%     VarianceMatrix = ones(size(RHSMatrix));
%     % If no variance matrix is not externally provided; consider equal
%     % weight for all points
% end
VarianceMatrix = ones(size(VarianceMatrix));

LHSMatrix = delY_delt_GPMat;
RHSMatrix = eqnParam*del2Y_delz2_GPMat;

% Trimming ~~~~
fraction_of_Zdata_removed = DataTrim.z; % Value in fraction in one Side
fraction_of_tdata_removed = DataTrim.t; % Value in fraction in one Side

% Twice of fraction_of_Zdata_removed will be removed

[ZZZ, TTT] = size(delY_delt_GPMat);
zdataselection = (ceil(fraction_of_Zdata_removed*ZZZ)+1):...
                (floor((1-fraction_of_Zdata_removed)*ZZZ));

tdataselection = (ceil(fraction_of_tdata_removed*TTT)+1):...
                (floor((1-fraction_of_tdata_removed)*TTT));           

% ############## Defining Error ####################
errMatrix = LHSMatrix - RHSMatrix;

errNorm1 = (sum(sum(abs(errMatrix))));
errNorm2 = sqrt(sum(sum((errMatrix(zdataselection, tdataselection).^2)./...
                                VarianceMatrix(zdataselection, tdataselection))));
errNormInf = max(max(abs(errMatrix)));