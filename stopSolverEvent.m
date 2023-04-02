function [value,stop,params] = stopSolverEvent(Y_all,stopParams,params)
    % The function that stops the simulation. The event functions doesn't
    % let us use the previouse step's Ys and only the Y for the current
    % step is usable in an event function. To use the previouse step's Ys
    % we had to use global variables but in the new update, we solved this
    % problem by not using event functions! we stop the simulation directly
    % from ode15s. This function is used inside ode15s_modified function
    % inside the solver folder.

    if 2 < size(Y_all,1)
        if ~ params.stop_th1
            [params.peaks_th1, params.peakCounter_th1, params.peakClass_th1, params.classCounter_th1, params.steadyStateAmp_th1, params.peaksInbetween_th1, params.lastCycleTime_th1, params.stop_th1] = ...
                stopSim(Y_all(end, 1),Y_all,size(Y_all,1)-1,stopParams.targetIndex_th1,stopParams.peakRepeat_th1,params.peakClass_th1,stopParams.peakTollerance,params.classCounter_th1,params.peakCounter_th1,params.peaks_th1,params.peaksInbetween_th1);
        end

        if ~ params.stop_th2
            [params.peaks_th2, params.peakCounter_th2, params.peakClass_th2, params.classCounter_th2, params.steadyStateAmp_th2, params.peaksInbetween_th2, params.lastCycleTime_th2, params.stop_th2] = ...
                stopSim(Y_all(end, 1),Y_all,size(Y_all,1)-1,stopParams.targetIndex_th2,stopParams.peakRepeat_th2,params.peakClass_th2,stopParams.peakTollerance,params.classCounter_th2,params.peakCounter_th2,params.peaks_th2,params.peaksInbetween_th2);
        end
    end

    if params.stop_th1 && params.stop_th2 
        % We dont want absorber's oscillations determine when we have
        % reached the steady state, because only one variable reaching ss
        % is enough to determine of we are in SS mode. But for debugging
        % and calculating the last cycle of absorber (Viscouse damper or 
        % SMA damper) we find its time response peaks. Also we have used
        % isfield so that this code could be generalized to systems which
        % don't have xd.
        if isfield(params, "stop_xd")
            % To get the last cycle of xd, we dont use the scheme we used for
            % th1 and th2; After the symulation is over, we provide xd (From 
            % end of the simulation to the beginning). We find the peaks and as
            % soon as we 
            for j=size(Y_all,1):-1:1

                [params.peaks_xd, params.peakCounter_xd, params.peakClass_xd, params.classCounter_xd, params.steadyStateAmp_xd, params.peaksInbetween_xd, params.lastCycleTime_xd, params.stop_xd] = ...
                    stopSim(Y_all(j-1, 1),Y_all,j-1,stopParams.targetIndex_xd,stopParams.peakRepeat_xd,params.peakClass_xd,stopParams.peakTollerance,params.classCounter_xd,params.peakCounter_xd,params.peaks_xd,params.peaksInbetween_xd);
                
                if params.stop_xd          
                    params.peaks_xd = params.peaks_xd(params.peaks_xd(:,1)~=0,:);
                    params.lastCycleTime_xd = params.peaks_xd(1:params.peaksInbetween_xd+1,1);
                    break;
                elseif j == 2
                    % No steady last cycle found
                    params.lastCycleTime_xd = [Y_all(end,1) Y_all(1,1)];
                    break;
                end
            end

            % Finding the upper and lower bound of the SS response. We
            % continue the search from end of the simulation untill the 
            % minimum of lastCycle of th1 or th2. We do so to avoid
            % processing the transient response
            upperBnd = params.peaks_xd(1,2);
            lowerBnd = params.peaks_xd(1,2);
            t_ = min(params.lastCycleTime_th1(1), params.lastCycleTime_th2(1));
            for i=1:length(params.peaks_xd)
                if upperBnd < params.peaks_xd(i,2)
                    upperBnd = params.peaks_xd(i,2);
                elseif params.peaks_xd(i,2) < lowerBnd
                    lowerBnd = params.peaks_xd(i,2);
                end

                % Break the loop.
                if params.peaks_xd(i,1) < t_
                    break;
                end
            end

            % Finding the pattern. The logic is as follows: First we have
            % to determine the bounds of the steady state response of xd.
            % After doing this we iterate the peaks of xd from end to start
            % and tag each peak as "U" or "L" for upper bound and lower
            % bound respectivels. The peaks that do not belong to upper
            % or lower bound, are assigned "-". The first peak that is not
            % tagged with "-" will be the firstBnd and the start of the last 
            % cycle. We will continue the search untill we reach a peak
            % that tagged with firstBnd's tag, But before that we should
            % have passed a peak with a tag, upposite of the firstBnd's
            % tag (By this, the variable triggle is changed from false to true).
            firstBnd = "";
            lastCycle = [Y_all(end,1) Y_all(1,1)];
            trigger = false;
            
            for i=1:length(params.peaks_xd)
                t_ = determineBound(params.peaks_xd(i,2), upperBnd, lowerBnd);
                if t_ ~= "-" && firstBnd == "" 
                    firstBnd = t_;
                    lastCycle(1) = params.peaks_xd(i,1);
                end
                
                if t_ ~= firstBnd && t_ ~= "-"
                    trigger = true;
                end

                if t_ == firstBnd && trigger
                    lastCycle(2) = params.peaks_xd(i,1);
                    params.lastCycleTime_xd = lastCycle;
                    break;
                end
%                 fprintf("%.4f || %.7f || %s\n", params.peaks_xd(i,1), params.peaks_xd(i,2), determineBound(params.peaks_xd(i,2), upperBnd, lowerBnd))
            end
        end


        % Set the start and end time of the last steadyState cycle
        if params.lastCycleTime_th1(end) < params.lastCycleTime_th2(end)
            r_ = params.lastCycleTime_th2;
        else
            r_ = params.lastCycleTime_th1;
        end

        if 2 <= size(r_,1)
            params.lastCycleTime = [r_(1), r_(end)];
        else
            params.lastCycleTime = [0 0];
        end

        % Set the values necessary for stopping the simulation
        value = 1;
        stop = 1;
    else
        params.lastCycleTime = [0 0];
        value = 0;
        stop = 0;
    end
    
    function str = determineBound(peak, upperBnd, lowerBnd)
        % This function determines if a given peak is for the upper bound
        % of the steady state response of for the lower bound. The
        % criterion of closeness is %2 of the benchmark peak.

        % Args:
        % peak: The amplitude of the peak
        % upperBnd: The upper bound of the steady state response
        % lowerBnd: The lower bound of the steady state response

        % Returns:
        % str: "U" if for upper bound, "L" for lower bund and "-" for a
        % pear that is not for upper nor lower bound.
        
        if abs(abs(peak) - abs(upperBnd))/abs(upperBnd) <= 0.02
            str = "U";
        elseif abs(abs(peak) - abs(lowerBnd))/abs(lowerBnd) <= 0.02
            str = "L";
        else
            str = "-";
        end
    end
end