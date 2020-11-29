function R_square = CoeffOfDeter(ObservedVWC, PredictedVWC)
    ObservedVWC_bar = nanmean(ObservedVWC);
    SStot = nansum((ObservedVWC-ObservedVWC_bar).^2);
    SSerr = nansum((ObservedVWC-PredictedVWC).^2);
    R_square = 1-(SSerr/SStot);
end