function RMSE = RootMeanSqrErr(ObservedVWC, PredictedVWC)
    Error = PredictedVWC - ObservedVWC;
    MeanSquareError = nanmean(Error.^2);
    RMSE = sqrt(MeanSquareError);
end