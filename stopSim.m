function [peaks, peakCounter, peakClass, classCounter, steadyStateAmp, peaksInbetween, lastCycleTime, stop] = stopSim(t,Y,i,targetIndex,peakRepeat,peakClass,peakTollerance,classCounter,peakCounter,peaks,peaksInbetween)
    % this function stops the simulation befor the simulation time ends. If
    % the harmonic response is periodic, there is no need to continue the
    % simulation because no new computations will occur. This function
    % looks for patternis in the response of the system and determines if
    % simulation can end.

    % Args:
    % t: simulation time
    % Y: a matrix consisting state variables in different times
    % (iterations) the first columns is time, the rest are the state
    % variables e.g. [time x dx]
    % i: the current line of the Y (Note that Y is a predefined matrix 
    % that updates as the simulation proceeds)
    % targetIndex: The number of column in Y matrix that we want to find
    % patterns in.
    %                        hoselam sar raft bagiyasho benevisam!!!!
    stop = false;
    steadyStateAmp = 0;
    lastCycleTime = 0;

    if sign(Y(i+1,targetIndex+1)) ~= sign(Y(i,targetIndex+1))
        peaks(peakCounter,:) = [t Y(i,targetIndex)];
        % Get the peaks that are close to the peak that we just found. If
        % the number of occurances are more than "peakrepeat - 1", add this
        % as a peak class to peakClass matrix. The first column of this
        % matrix is the peak's amplitude, the rest of the columns are the
        % peak number of the occurances.
        x = find( ...
            abs(peaks(1:peakCounter,2)-peaks(peakCounter,2))/abs(peaks(peakCounter,2))< peakTollerance, ...
            peakRepeat, "last");

        x = x(1:end-1);
        if(peakRepeat-1 <= length(x))
            
            xx = find( abs(peakClass(:,1) - Y(i,targetIndex))/abs(Y(i,targetIndex)) < peakTollerance, 1);
            if ~isempty(xx)
                if peaksInbetween == 0
                    % Get the number of peaks in each period
                    peaksInbetween = peakClass(classCounter-1,peakRepeat) - peakClass(classCounter-1,peakRepeat-1);
                end
                
                % Break the loop
                steadyStateAmp = max(abs(peaks(peakCounter-peaksInbetween:peakCounter,2)));
                lastCycleTime = peaks(peakCounter-peaksInbetween:peakCounter,1);

                stop = true;
            end

            peakClass(classCounter,:) = [Y(i,targetIndex) x'];
            classCounter = classCounter + 1;

        end

        peakCounter = peakCounter + 1;
    end

end
