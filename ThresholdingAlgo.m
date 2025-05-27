% This function is from: Brakel, J.P.G. van (2014). "Robust peak detection algorithm using z-scores". Stack Overflow. Available at: https://stackoverflow.com/questions/22583391/peak-signal-detection-in-realtime-timeseries-data/22640362#22640362 (version: 2020-11-08).

function [signals,avgFilter,stdFilter] = ThresholdingAlgo(y,lag,threshold,influence,T) %had, Cell also after ",influence"
%Initialise signal results
signals=zeros(length(y),1);
%signals = (length(y),Cell);
% Initialise filtered series
filteredY = y(1:lag+1);        %Take the n (lag+1) first value of the trace 
% Initialise filters
avgFilter(lag+1,1) = mean(y(1:lag+1));  %do the mean of the n(lag+1) first value of the trace
stdFilter(lag+1,1) = std(y(1:lag+1));   %do the std of the n(lag+1) first value of the trace
% Loop over all datapoints y(lag+2),...,y(t)
for i=lag+2:length(y)
    % If new value is a specified number of deviations away
    if abs(y(i)-avgFilter(i-1)) > threshold*stdFilter(i-1)   %If absolute value of (fluo value - avrg filter) is superior of std filtred times threshold
        if y(i) > avgFilter(i-1)         %If fluo value superior to avrg filter value
            % Positive signal
            signals(i) = 1;                                        % Then Signal =1 
        else
            % Negative signal 
            signals(i) = -1;                                       %Otherwise signal=-1
        end
        % Make influence lower
        filteredY(i) = influence*y(i)+(1-influence)*filteredY(i-1); 
    else
        % No signal
        signals(i) = 0;                                      %Else if fluo value - avrg isn't superior to thr*std, there's no signal and signal=0
        filteredY(i) = y(i);                                 % Put one more value of y in the filtered vector
    end
    % Adjust the filters
    avgFilter(i) = mean(filteredY(i-lag:i));   %Shift avrg Filter and adjust them with the influence if there's some
    stdFilter(i) = std(filteredY(i-lag:i));    %Shift std Filter
end  
% Done, now return results
end





