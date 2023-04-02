function [frequencyRange, amplitude, savedData] = HandSystem_3DOF_AbsorberOnElbow_CODE_MULTIPLE_PARALLEL(verbose, plotFRF, saveData, sysObj, options)
    nOutVars = nargout; % The number of variables requested to return
    % options: struct: should have errToll and simStep. If it has a
    % message, at the beginning of each iteration, it will be displayed
    
    % savedData: struct: It is a struct that will keep the SMA spring's
    % tension-displacement data. the data for each frequency can be
    % accessed by "iter_<frequency>" field name. It will contain the necessary
    % displacements and force/stress. The choice to save stress of force of
    % the spring in each step is determined in ode15s file. To change it,
    % refer to the "sovler" folder. Also if no SMA spring is used in the
    % system, it will stay empty.
    
    % Start a parallel pool if non exists, make one
    if verbose
        pool = gcp;
        fprintf("Allocated %d cores.\n", pool.NumWorkers);
    end
    
    system = sysObj;
    % The element that will change on every iteration
    freqSweepe = system.freqSweepe;

    saveBucket = zeros(length(freqSweepe), 3);
    savedData_ = cell(length(freqSweepe),1); % The variable that contains the SMA data. (Null if no spring used)

    % Locating the solver folder and adding it to path
    addpath(sprintf("./solver/R2019b"))
    parfor iter = 1:length(freqSweepe)
        odeParams = struct(); % The variable that will be passed into ode15s
        dynamicSysParams = struct(); % The variable that will be passed into dynamic model
        stopParams = struct(); % The variable that will be passed into  ode15s but will be used in stopSim function

        dynamicSysParams.m1 = system.m1;
        dynamicSysParams.m2 = system.m2;
        dynamicSysParams.md = system.md;
        dynamicSysParams.g = system.g;
        dynamicSysParams.c1 = system.c1;
        dynamicSysParams.c2 = system.c2;
        dynamicSysParams.kt1Mul = system.kt1Mul;
        dynamicSysParams.kt2Mul = system.kt2Mul;
        dynamicSysParams.l1 = system.l1;
        dynamicSysParams.l2 = system.l2;
        dynamicSysParams.lc1 = system.lc1;
        dynamicSysParams.lc2 = system.lc2;
        dynamicSysParams.ld = system.ld;
        dynamicSysParams.xdInit = system.xdInit;
        dynamicSysParams.i1 = system.i1;
        dynamicSysParams.i2 = system.i2;
        dynamicSysParams.cd = system.cd;
        dynamicSysParams.kd = system.kd;
        dynamicSysParams.spring1_Active = system.spring1.Active;
        dynamicSysParams.spring2_Active = system.spring2.Active;
    
        dynamicSysParams.ampWrist = system.ampWrist;
        dynamicSysParams.ampElbow = system.ampElbow;

        steadyStateAmp_th1 = 0;
        steadyStateAmp_th2 = 0;
        steadyStateAmp_xd = 0;
    
        % The initial state of the system
        th1_Init = 0;
        Dth1_Init = 0;
        th2_Init = 0;
        Dth2_Init = 0;
        xd_Init = 0;
        Dxd_Init = 0;
        
        Y = [th2_Init,Dth2_Init,th1_Init,Dth1_Init,xd_Init,Dxd_Init];
        
        % ODE integration ime bound and increments
        tStart = 0;
        tEnd = 120;
        increment = options.simStep;


        if verbose; clc; fprintf("%s\n", options.msg); end
        if verbose; fprintf("Freq: %.2f\n", freqSweepe(iter)); end

        % If frequency of a DOF is -1, it means that we want it to use
        % frequency sweepe property as its value throughout the simulation
        if system.freqWrist == -1; dynamicSysParams.freqWrist = freqSweepe(iter); else; dynamicSysParams.freqWrist = 0; end
        if system.freqElbow == -1; dynamicSysParams.freqElbow = freqSweepe(iter); else; dynamicSysParams.freqElbow = 0; end


        % defining the parameters necessary for spring1we instantiate these
        % variables so that we can use ode15s (modified version) 
        if dynamicSysParams.spring1_Active
            spring1 = struct(); spring1_temp = struct();
            z0 = 0; zS0 = 0; zT0 = 0; z = 0;
            eDot = 1; phase = 1; prevPhase = 1; n1 = 10000;
            dynamicSysParams.z0 = z0; dynamicSysParams.zS0 = zS0; dynamicSysParams.forceLinear = false; dynamicSysParams.forceLinearZ0 = nan; dynamicSysParams.forceLinearStress = nan; dynamicSysParams.forceLinearDir = nan;
            spring1_temp = SpringClass().init(system.spring1.d0,system.spring1.D0,system.spring1.l0,system.spring1.N,system.spring1.S_M_finish,system.spring1.S_M_start,system.spring1.S_A_start,system.spring1.S_A_finish,system.spring1.Gm,system.spring1.Ga,system.spring1.eL,system.spring1.T,system.spring1.As,system.spring1.Af,system.spring1.Ms,system.spring1.Mf,z0,zS0,zT0,z,eDot,phase,prevPhase,n1,system.spring1.preStretchPercentage,dynamicSysParams, system.spring1.name);
            spring1_temp.verbose = false;
        
            dynamicSysParams.spring1 = spring1;
            dynamicSysParams.spring1_temp = spring1_temp;
        end
        
        % defining the parameters necessary for spring2we instantiate these
        % variables so that we can use ode15s (modified version) 
        if dynamicSysParams.spring2_Active
            spring2 = struct(); spring2_temp = struct();
            z0 = 0; zS0 = 0; zT0 = 0; z = 0;
            eDot = 1; phase = 1; prevPhase = 1; n1 = 10000;
            dynamicSysParams.z0 = z0; dynamicSysParams.zS0 = zS0; dynamicSysParams.forceLinear = false; dynamicSysParams.forceLinearZ0 = nan; dynamicSysParams.forceLinearStress = nan; dynamicSysParams.forceLinearDir = nan;
            spring2_temp = SpringClass().init(system.spring2.d0,system.spring2.D0,system.spring2.l0,system.spring2.N,system.spring2.S_M_finish,system.spring2.S_M_start,system.spring2.S_A_start,system.spring2.S_A_finish,system.spring2.Gm,system.spring2.Ga,system.spring2.eL,system.spring2.T,system.spring2.As,system.spring2.Af,system.spring2.Ms,system.spring2.Mf,z0,zS0,zT0,z,eDot,phase,prevPhase,n1,system.spring2.preStretchPercentage,dynamicSysParams, system.spring2.name);
            spring2_temp.verbose = false;
        
            dynamicSysParams.spring2 = spring2;
            dynamicSysParams.spring2_temp = spring2_temp;
        end


        % Defining the parameters necessary form event fucntion to stop the
        % integration
        peakRepeat_th1 = 10;
        peakRepeat_th2 = 10;
        peakRepeat_xd = 10;
        stopParams.peakDelay = 1; % After the stopping criterion is met, we delay the stopping of the system to be more sure of the steadiness of the amplitude. 
        stopParams.peakTollerance = options.peakErrToll; % The amount of peak reduction relative to its previous peak to stop the symulation
        stopParams.targetIndex_th1 = 3 + 1;
        stopParams.targetIndex_th2 = 1 + 1; 
        stopParams.targetIndex_xd = 5 + 1; 
        stopParams.peakRepeat_th1 = peakRepeat_th1;
        stopParams.peakRepeat_th2 = peakRepeat_th2;
        stopParams.peakRepeat_xd = peakRepeat_xd;
        
        odeParams.peakClass_th1 = zeros(100,peakRepeat_th1);
        odeParams.peakClass_th2 = zeros(100,peakRepeat_th2);
        odeParams.peakClass_xd = zeros(100,peakRepeat_xd);
        odeParams.peaks_th1 = zeros(200,2); odeParams.peakCounter_th1 = 1;
        odeParams.peaks_th2 = zeros(200,2); odeParams.peakCounter_th2 = 1;
        odeParams.peaks_xd = zeros(200,2); odeParams.peakCounter_xd = 1;
        odeParams.stop_th1 = false;
        odeParams.stop_th2 = false;
        odeParams.stop_xd = false;
        odeParams.classCounter_th1 = 1;
        odeParams.classCounter_th2 = 1;  
        odeParams.classCounter_xd = 1;    
        odeParams.steadyStateAmp_th1 = 0;
        odeParams.steadyStateAmp_th2 = 0;
        odeParams.steadyStateAmp_xd = 0;
        odeParams.peaksInbetween_th1 = 0;
        odeParams.peaksInbetween_th2 = 0;
        odeParams.peaksInbetween_xd = 0;
        odeParams.stopParams = stopParams;

        % ODE solver options
        odeOptions = odeset('RelTol', options.RelTol, 'AbsTol', options.AbsTol);

        try
            [odeParams,~,Y] = ode15s_Modified( ...
                @(t,Y,params) dynamicEqn_3DOF_AbsorberOnElbow(t,Y,params), ...
                linspace(tStart,tEnd,(tEnd-tStart)/increment), ...
                [th2_Init;Dth2_Init;th1_Init;Dth1_Init;xd_Init;Dxd_Init], odeOptions, dynamicSysParams, odeParams, dynamicSysParams);
            steadyStateAmp_th1 = odeParams.steadyStateAmp_th1;
            steadyStateAmp_th2 = odeParams.steadyStateAmp_th2;
            steadyStateAmp_xd = 0;

            if 1 < size(odeParams.savedData,1) && 2 < nOutVars
                % Only add to saveData if we have saved any datas while solving
                % the ode, else add nan. Also if savedData is not requested
                % at the retun, no need to fill the structure

                if dynamicSysParams.freqElbow ~= 0
                    % Get the lastCycle's data. The data is in the following
                    % format: [time, displacement, savedData1, savedData2]
                    tmp_ = [odeParams.accepted_local(:,1), odeParams.accepted_local(:,6), odeParams.savedData];
                    tmp_ = tmp_(odeParams.lastCycleTime_xd(end)<tmp_(:,1) & tmp_(:,1)<odeParams.lastCycleTime_xd(1),:);

                    % Getting the kinetic energy of the last cycle
                    odeParams_LastCycle = odeParams.accepted_local( ...
                        odeParams.lastCycleTime_xd(end)<odeParams.accepted_local(:,1) & odeParams.accepted_local(:,1) <odeParams.lastCycleTime_xd(1) ...
                        ,:);
                    th2 = odeParams_LastCycle(:,2); Dth2 = odeParams_LastCycle(:,3);
                    th1 = odeParams_LastCycle(:,4); Dth1 = odeParams_LastCycle(:,5);
                    xd = odeParams_LastCycle(:,6); Dxd = odeParams_LastCycle(:,7);
    
                    K_ =  ...
                        1/2*dynamicSysParams.m1*( (dynamicSysParams.lc1*Dth1).^2 ) + ...
                        1/2*dynamicSysParams.i1*(Dth1).^2 + ...
                        1/2*dynamicSysParams.i2*(Dth1+Dth2).^2 + ...
                        1/2*dynamicSysParams.md*( ...
                            (sqrt(xd.^2+dynamicSysParams.ld^2).*Dth1).^2 + ...
                            Dxd.^2 + ...
                            2*Dth1.*Dxd*dynamicSysParams.ld ...
                        ) + ...
                        1/2*dynamicSysParams.m2*( ...
                            (dynamicSysParams.l1*Dth1).^2 + ...
                            (dynamicSysParams.lc2*(Dth2+Dth1)).^2 + ... 
                            2*dynamicSysParams.l1*dynamicSysParams.lc2*Dth1.*(Dth2+Dth1).*cos(th2) ...
                        );
                    KineticEnergy_LastCycle = mean(K_);
                else
                    % If the excitation frequency is zero, There is no last
                    % cycle
                    tmp_ = [0 0 0 0];
                    KineticEnergy_LastCycle = 0;
                end

                % Add the data
                savedData_{iter} = struct( ...
                    "frequency", freqSweepe(iter), ...
                    "spring1", odeParams.savedData(:,1), ...
                    "spring2", odeParams.savedData(:,2),...
                    "displacement", odeParams.accepted_local(:,6), ...
                    "Ed", polyarea(tmp_(:,2),tmp_(:,3)) + polyarea(tmp_(:,2),tmp_(:,4)), ...
                    "lastCycle", tmp_, ...
                    "K_Last", KineticEnergy_LastCycle ...
                );
            elseif 2 < nOutVars
                % If we have no SMAs, calculate the energy dissipated in
                % each cycle through cd. First we find the last cycle of
                % xd. The logic is that we calculate integral(F*dxd) in the
                % last cycle. (xd is the displacement of the absorber and 
                % F = cd * Dxd)
                if dynamicSysParams.freqElbow ~= 0
                    % Add t, xd and Force excerted by damper to a matrix
                    % and slice the data to the last cycle.
                    tmp_ = [odeParams.accepted_local(:,1), odeParams.accepted_local(:,6), sysObj.cd*odeParams.accepted_local(:,7)];
                    tmp_ = tmp_(odeParams.lastCycleTime_xd(end)<tmp_(:,1) & tmp_(:,1)<odeParams.lastCycleTime_xd(1),:);
    
                    % Here we connect end of the cure to its start and calculate its area
                    Ed = polyarea([tmp_(:,2) ;tmp_(1,2)], [tmp_(:,3) ;tmp_(1,3)]); 
                    % Getting the kinetic energy of the last cycle
                    odeParams_LastCycle = odeParams.accepted_local( ...
                        odeParams.lastCycleTime_xd(end)<odeParams.accepted_local(:,1) & odeParams.accepted_local(:,1) <odeParams.lastCycleTime_xd(1) ...
                        ,:);
                    th2 = odeParams_LastCycle(:,2); Dth2 = odeParams_LastCycle(:,3);
                    th1 = odeParams_LastCycle(:,4); Dth1 = odeParams_LastCycle(:,5);
                    xd = odeParams_LastCycle(:,6); Dxd = odeParams_LastCycle(:,7);
    
                    K_ =  ...
                        1/2*dynamicSysParams.m1*( (dynamicSysParams.lc1*Dth1).^2 ) + ...
                        1/2*dynamicSysParams.i1*(Dth1).^2 + ...
                        1/2*dynamicSysParams.i2*(Dth1+Dth2).^2 + ...
                        1/2*dynamicSysParams.md*( ...
                            (sqrt(xd.^2+dynamicSysParams.ld^2).*Dth1).^2 + ...
                            Dxd.^2 + ...
                            2*Dth1.*Dxd*dynamicSysParams.ld ...
                        ) + ...
                        1/2*dynamicSysParams.m2*( ...
                            (dynamicSysParams.l1*Dth1).^2 + ...
                            (dynamicSysParams.lc2*(Dth2+Dth1)).^2 + ... 
                            2*dynamicSysParams.l1*dynamicSysParams.lc2*Dth1.*(Dth2+Dth1).*cos(th2) ...
                        );
                    KineticEnergy_LastCycle = mean(K_);
                else 
                    % If the excitation frequency is zero, we assume no
                    % energy absorbtion.
                    Ed = 0;
                    KineticEnergy_LastCycle = 0;
                end

                % Add the data
                savedData_{iter} = struct( ...
                    "frequency", freqSweepe(iter), ...
                    "Ed", Ed, ...
                    "K_Last", KineticEnergy_LastCycle ...
                );
            else
                savedData_{iter} = nan;
            end
        catch msg
            disp(msg.message)
            steadyStateAmp_th1 = 999;
            steadyStateAmp_th2 = 999;
            steadyStateAmp_xd  = 999;
            
            % Add nan if we ran to errors.
            savedData_{iter} = nan;
        end
        
        if steadyStateAmp_th1 == 0
            % If there are no steady state amps found, save the maximum of the last
            % 10 peaks as steady state amplitude
            peaks_th1 = odeParams.peaks_th1;
            peaks_th1 = peaks_th1(peaks_th1(:,2) ~= 0,:);
            if verbose; disp("No steady statate amplitudes were found for th1"); end
        
            if 0 < length(peaks_th1)-10
                steadyStateAmp_th1 = max(abs(  peaks_th1(length(peaks_th1)-10:end,2)  ));
            else
                % If we have less than 10 peaks
                steadyStateAmp_th1 = max(abs(  peaks_th1(:,2)  ));
            end
        end
        
        if steadyStateAmp_th2 == 0
            % If there are no steady state amps found, save the maximum of the last
            % 10 peaks as steady state amplitude
            peaks_th2 = odeParams.peaks_th2;
            peaks_th2 = peaks_th2(peaks_th2(:,2) ~= 0,:);
            if verbose; disp("No steady statate amplitudes were found for th2"); end
        
            if 0 < length(peaks_th2)-10
                steadyStateAmp_th2 = max(abs(  peaks_th2(length(peaks_th2)-10:end,2)  ));
            else
                % If we have less than 10 peaks
                steadyStateAmp_th2 = max(abs(  peaks_th2(:,2)  ));
            end
        end

        % finding the amplitude of absorber
        out = Y(:,5);
        xdPeaks = [];
        for i = length(out)-1:-1:0
            if 1<i
                if out(i-1) < out(i) && out(i+1) < out(i)
                    xdPeaks(end+1) = out(i); 
                elseif out(i-1) > out(i) && out(i+1) > out(i)
                    xdPeaks(end+1) = out(i);
                end
            
                if length(xdPeaks) >= 2
                    steadyStateAmp_xd = 1000*max(abs(xdPeaks(1)) , abs(xdPeaks(2)));
                    break;
                end
            else
                steadyStateAmp_xd = 0;
            end
        end
        
        saveBucket(iter,:) = [steadyStateAmp_th1 steadyStateAmp_th2 steadyStateAmp_xd];
    end
    savedData = savedData_;

    % Removing the solver from matlab's path
    rmpath(sprintf("./solver/R2019b"))

    if plotFRF
        hold on; legend show;
        plotFRF(freqSweepe, saveBucket, "DisplayName", strcat("Damped:: ld:", num2str(ld), " | kd:", num2str(kd), " | cd:", num2str(cd) ));
    end

    frequencyRange = freqSweepe';

    % Saving the amplitude daya as a structure for better readability
    amplitude = struct();
    amplitude.th1 = saveBucket(:,1);
    amplitude.th2 = saveBucket(:,2);
    amplitude.xd = saveBucket(:,2);

    if saveData
        writematrix([frequencyRange amplitude], sprintf("./Datas/Damped_ld_%.5f#kd_%.5f#cd_%.5f#.csv", ld, kd, cd));
    end

end
