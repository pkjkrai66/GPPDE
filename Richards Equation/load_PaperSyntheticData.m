function [zVector, tVector, SM_Observation, thta, Ks, alpha, N, hyper0 ] = ...
                                      load_PaperSyntheticData(DataSetName)
                                    
                                    
% Options for DataSetName:
%      SyntheticSMData_1ErrorFree
%      SyntheticSMData_1ErrorFree_12Points
%      SyntheticSMData_1Errored
%      SyntheticSMData_1Errored_12Points
switch DataSetName
    case 'SyntheticSMData_1ErrorFree'
        load('SyntheticSMData_1ErrorFree.mat')
        hyper0 = [5, 0.5, 0.001];       % Converges to -->
        % --> [3.42831211976612, 0.624365255593804, 0.0000128384776762633]

    case 'SyntheticSMData_1ErrorFree_12Points'
        load('SyntheticSMData_1ErrorFree.mat')
        zVector = z.discreteZ([1, [6 12 17 23 28 34 39 45 50 56], end]);
        SM_Observation2 = SM_Observation;
        SM_Observation = SM_Observation2(...
                            [1, [6 12 17 23 28 34 39 45 50 56], end],:);
        hyper0 = [10 50 0.02];       % Converges to -->
        % --> [10.0165958651257, 49.9930563158343, -0.000448685635114246]
        
    case 'SyntheticSMData_1Errored'
        load('SyntheticSMData_1Errored')
        hyper0 = [6, 0.02, 0.002];      % Converges to -->
        % --> [5.99040717875626, 0.0264296285299399, 0.0025845719840390]

    case 'SyntheticSMData_1Errored_12Points'
        load('SyntheticSMData_1Errored')
        zVector = z.discreteZ([1, [6 12 17 23 28 34 39 45 50 56], end]);
        SM_Observation2 = SM_Observation;
        SM_Observation = SM_Observation2(...
                        [1, [6 12 17 23 28 34 39 45 50 56], end],:);
        hyper0 = [6, 0.02, 0.002];      % Converges to -->
        % --> [7.48275974340419,0.0261691154321685,0.00267869348221790]
end
