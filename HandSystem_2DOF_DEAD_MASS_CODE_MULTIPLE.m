function [frequencyRange, amplitude, savedData] = HandSystem_2DOF_DEAD_MASS_CODE_MULTIPLE(verbose, plotFRF, saveData, sysObj, options)
    % options: struct: should have errToll and simStep. If it has a
    % message, at the beginning of each iteration, it will be displayed

    odeParams = struct(); % The variable that will be passed into ode15s
    dynamicSysParams = struct(); % The variable that will be passed into dynamic model
    stopParams = struct(); % The variable that will be passed into  ode15s but will be used in stopSim function

    system = sysObj;
    % The element that will change on every iteration
    freqSweepe =  system.freqSweepe;

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

    dynamicSysParams.ampWrist = system.ampWrist;
    dynamicSysParams.ampElbow = system.ampElbow;

    % Necessary for not running into errors
    dynamicSysParams.spring1_Active = false;
    dynamicSysParams.spring2_Active = false;

    % The initial state of the system
    th1_Init = 0;
    Dth1_Init = 0;
    th2_Init = 0;
    Dth2_Init = 0;

    saveBucket = zeros(length(freqSweepe), 2);
    savedData_ = cell(length(freqSweepe),1); % The variable that contains the SMA data. (Null if no spring used)

    % ODE integration ime bound and increments
    tStart = 0;
    tEnd = 120;
    increment = options.simStep;
    
    % Locating the solver folder and adding it to path
    addpath(sprintf("./solver/R%s/", "2019b"))
    for iter = 1:1:length(freqSweepe)

        if verbose; clc; fprintf("%s \n", options.msg); end
        if verbose; fprintf("Freq: %.1f\n", freqSweepe(iter)); end

        % If frequency of a DOF is -1, it means that we want it to use
        % frequency sweepe property as its value throughout the simulation
        if system.freqWrist == -1; dynamicSysParams.freqWrist = freqSweepe(iter); else; dynamicSysParams.freqWrist = 0; end
        if system.freqElbow == -1; dynamicSysParams.freqElbow = freqSweepe(iter); else; dynamicSysParams.freqElbow = 0; end
        
        % Defining the parameters necessary form event fucntion to stop the
        % integration
        peakRepeat_th1 = 4;
        peakRepeat_th2 = 4;
        stopParams.peakDelay = 1; % After the stopping criterion is met, we delay the stopping of the system to be more sure of the steadiness of the amplitude. 
        stopParams.peakTollerance = options.peakErrToll; % The amount of peak reduction relative to its previous peak to stop the symulation
        stopParams.targetIndex_th1 = 3 + 1;
        stopParams.targetIndex_th2 = 1 + 1; 
        stopParams.peakRepeat_th1 = 4;
        stopParams.peakRepeat_th2 = 4;
        
        odeParams.peakClass_th1 = zeros(100,peakRepeat_th1);
        odeParams.peakClass_th2 = zeros(100,peakRepeat_th2);
        odeParams.peaks_th1 = zeros(200,2); odeParams.peakCounter_th1 = 1;
        odeParams.peaks_th2 = zeros(200,2); odeParams.peakCounter_th2 = 1;
        odeParams.stop_th1 = false;
        odeParams.stop_th2 = false;
        odeParams.classCounter_th1 = 1;
        odeParams.classCounter_th2 = 1;    
        odeParams.steadyStateAmp_th1 = 0;
        odeParams.steadyStateAmp_th2 = 0;
        odeParams.peaksInbetween_th1 = 0;
        odeParams.peaksInbetween_th2 = 0;
        odeParams.stopParams = stopParams;

        % ODE solver options
        odeOptions = odeset('RelTol', options.RelTol, 'AbsTol', options.AbsTol);
        
        [odeParams,~,~] = ode15s_Modified( ...
            @(t,Y,params) dynamicEqn_2DOF_DEAD_MASS(t,Y,params), ...
            linspace(tStart,tEnd,(tEnd-tStart)/increment), ...
            [th2_Init;Dth2_Init;th1_Init;Dth1_Init], odeOptions, dynamicSysParams, odeParams, dynamicSysParams);

        steadyStateAmp_th1 = odeParams.steadyStateAmp_th1;
        steadyStateAmp_th2 = odeParams.steadyStateAmp_th2;
        
        if dynamicSysParams.freqElbow ~= 0
            % Get the lastCycle's data. The data is in the following
            % format: [time, displacement, savedData1, savedData2]
            tmp_ = [odeParams.accepted_local(:,1), odeParams.accepted_local(:,2)];
            tmp_ = tmp_(odeParams.lastCycleTime_th2(1)<tmp_(:,1) & tmp_(:,1)<odeParams.lastCycleTime_th2(end),:);

            % Getting the kinetic energy of the last cycle
            odeParams_LastCycle = odeParams.accepted_local( ...
                odeParams.lastCycleTime_th1(1)<odeParams.accepted_local(:,1) & odeParams.accepted_local(:,1) <odeParams.lastCycleTime_th1(end) ...
                ,:);
            th2 = odeParams_LastCycle(:,2); Dth2 = odeParams_LastCycle(:,3);
            th1 = odeParams_LastCycle(:,4); Dth1 = odeParams_LastCycle(:,5);

            K_ =  ...
                1/2*dynamicSysParams.m1*( (dynamicSysParams.lc1*Dth1).^2 ) + ...
                1/2*dynamicSysParams.i1*(Dth1).^2 + ...
                1/2*dynamicSysParams.i2*(Dth1+Dth2).^2 + ...
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
            "K_Last", KineticEnergy_LastCycle ...
        );



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
        
        saveBucket(iter,:) = [steadyStateAmp_th1 steadyStateAmp_th2];
    end
    savedData = savedData_;
    % Removing the solver from matlab's path
    rmpath(sprintf("./solver/R%s/", "2019b"))

    if plotFRF
        hold on; legend show;
        plotFRF(freqSweepe, saveBucket, "DisplayName", strcat("Damped:: ld:", num2str(ld)));
    end

    frequencyRange = freqSweepe';

    % Saving the amplitude daya as a structure for better readability
    amplitude.th1 = saveBucket(:,1);
    amplitude.th2 = saveBucket(:,2);

    if saveData
        writematrix([frequencyRange amplitude], sprintf("./Datas/Damped_ld_%.5f#kd_%.5f#cd_%.5f#.csv", ld, kd, cd));
    end

end
