clc
clear all;
system = sysParams();
system = system.init();

odeParams = struct(); % The variable that will be passed into ode15s
dynamicSysParams = struct(); % The variable that will be passed into dynamic model
stopParams = struct(); % The variable that will be passed into  ode15s but will be used in stopSim function

dynamicSysParams.m1 = system.m1;
dynamicSysParams.m2 = system.m2;
dynamicSysParams.g  =  system.g;
dynamicSysParams.md = system.md;
dynamicSysParams.c1 = system.c1;
dynamicSysParams.c2 = system.c2;
dynamicSysParams.kt1Mul = system.kt1Mul;
dynamicSysParams.kt2Mul = system.kt2Mul;
dynamicSysParams.cd = system.cd;
dynamicSysParams.kd = system.kd;
dynamicSysParams.l1 = system.l1;
dynamicSysParams.l2 = system.l2;
dynamicSysParams.lc1 = system.lc1;
dynamicSysParams.lc2 = system.lc2;
dynamicSysParams.ld = system.ld;
dynamicSysParams.xdInit = system.xdInit;
dynamicSysParams.i1 = system.i1;
dynamicSysParams.i2 = system.i2;

dynamicSysParams.ampElbow = system.ampElbow;
dynamicSysParams.freqElbow = system.freqElbow;
dynamicSysParams.ampWrist = system.ampWrist;
dynamicSysParams.freqWrist= system.freqWrist;

dynamicSysParams.spring1_Active = system.spring1.Active;
dynamicSysParams.spring2_Active = system.spring2.Active;

%%%%%%%%%%%%%%%%%%
%%% Override system parameters

dynamicSysParams.ampElbow = 1;
dynamicSysParams.ampWrist = 0;
dynamicSysParams.freqElbow = 10;
dynamicSysParams.freqWrist = 10;

dynamicSysParams.md = .15;
dynamicSysParams.kd = 100;          
dynamicSysParams.cd = 1 ;       
dynamicSysParams.ld = 0.2;

dynamicSysParams.c1 = .1;
dynamicSysParams.c2 = .02;
dynamicSysParams.kt1Mul = 1;1/2.9564*6;
dynamicSysParams.kt2Mul = 1;1/1.6781*3;
 
dynamicSysParams.spring1_Active = false;
dynamicSysParams.spring2_Active = false;
%%%%%%%%%%%%%%%%%%

% The initial state of the system
th1_Init = 0;
Dth1_Init = 0;
th2_Init = 0;
Dth2_Init = 0;
xd_Init = 0;
Dxd_Init = 0;

% % For debugging purposes
% global accepted_local_G savedData_G

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
peakRepeat_th1 = 4;
peakRepeat_th2 = 4;
peakRepeat_xd = 4;
stopParams.peakDelay = 1; % After the stopping criterion is met, we delay the stopping of the system to be more sure of the steadiness of the amplitude. 
stopParams.peakTollerance = 0.005; % The amount of peak reduction relative to its previous peak to stop the symulation
stopParams.targetIndex_th1 = 3 + 1;
stopParams.targetIndex_th2 = 1 + 1; 
stopParams.targetIndex_xd = 5 + 1; 
stopParams.peakRepeat_th1 = 4;
stopParams.peakRepeat_th2 = 4;
stopParams.peakRepeat_xd = 4;

odeParams.peakClass_th1 = zeros(100,peakRepeat_th1);
odeParams.peakClass_th2 = zeros(100,peakRepeat_th2);
odeParams.peakClass_xd = zeros(100,peakRepeat_th2);
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

tStart = 0;
tEnd = 120;
increment = 0.01;

odeOptions = odeset('RelTol', 1e-5, 'AbsTol', 1e-7);

% Locating the solver folder and adding it to path
addpath(sprintf("./solver/R%s/", "2019b"))
tic
% Note that ode15s_Modified is teh mdoified version of ode15s. Also if you
% notice we have assigned dynamicSysParams twice for the function; the last
% one is a varargin and is imperative to updating the SMA spring state
% variables.
[odeParams,t,Y] = ode15s_Modified( ...
    @(t,Y,params) dynamicEqn_3DOF_AbsorberOnElbow(t,Y,params), ...
    linspace(tStart,tEnd,(tEnd-tStart)/increment), ...
    [th2_Init;Dth2_Init;th1_Init;Dth1_Init;xd_Init;Dxd_Init], odeOptions, dynamicSysParams, odeParams, dynamicSysParams);
toc
% Removing the solver from matlab's path
rmpath(sprintf("./solver/R%s/", "2019b"))

steadyStateAmp_th1 = odeParams.steadyStateAmp_th1;
steadyStateAmp_th2 = odeParams.steadyStateAmp_th2;

if steadyStateAmp_th1 == 0
    % If there are no steady state amps found, save the maximum of the last
    % 10 peaks as steady state amplitude
    peaks_th1 = odeParams.peaks_th1;
    peaks_th1 = peaks_th1(peaks_th1(:,2) ~= 0,:);
    disp("No steady statate amplitudes were found for th1")

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
    disp("No steady statate amplitudes were found for th2")

    if 0 < length(peaks_th2)-10
        steadyStateAmp_th2 = max(abs(  peaks_th2(length(peaks_th2)-10:end,2)  ));
    else
        % If we have less than 10 peaks
        steadyStateAmp_th2 = max(abs(  peaks_th2(:,2)  ));
    end
end

out = [];
out.tout = t; 
out.th2 = Y(:,1); 
out.Dth2= Y(:,2); 
out.th1 = Y(:,3); 
out.Dth1 = Y(:,4); 
out.xd = Y(:,5); 
out.Dxd = Y(:,6);
out.absorberForce = out.xd*dynamicSysParams.kd + out.Dxd*dynamicSysParams.cd;

f = figure; hold on; grid off; plot(out.tout, out.th1); f.Position = [100 100 600 400];
if dynamicSysParams.spring1_Active || dynamicSysParams.spring2_Active
    onlyLastCycle = true;

    if onlyLastCycle
        tmp = [odeParams.accepted_local(:,1), odeParams.accepted_local(:,6), odeParams.savedData];
        tmp = tmp(odeParams.lastCycleTime(1)<tmp(:,1) & tmp(:,1)<odeParams.lastCycleTime(2),:);
        f = figure; 
        hold on; grid off; xlabel("xd (m)"); ylabel("Force (N)"); xline(0, "--",'HandleVisibility','off'); yline(0, "--",'HandleVisibility','off');
        plot(tmp(:,2), tmp(:,3), "DisplayName", "spring 1");  f.Position = [800 100 400 400];
        if dynamicSysParams.spring2_Active
            hold on; plot(-tmp(:,2), tmp(:,4), "DisplayName", "spring 2");
        end
        legend show;
        legend('Location','northwest')
    else
        SS_SavedData = odeParams.savedData;
        SS_accepted = odeParams.accepted_local;
    
        % Deleting the data for transient answer
        limit = 90000000000;
        delStart = find((limit<SS_SavedData(:,2)),1,"first");
        delEnd = find((limit<SS_SavedData(:,2)),1,"last");
    
        % Trims a specific % of the beginning of the data
        trimData = 0; trimData = trimData / 100;
        if ~isempty(delEnd) && ~isempty(delStart) && trimData == 0
            SS_SavedData = [SS_SavedData(1:delStart,:); SS_SavedData(delEnd:end,:)];
            SS_accepted = [SS_accepted(1:delStart,:); SS_accepted(delEnd:end,:)];
        elseif trimData ~= 0
            SS_SavedData = SS_SavedData(ceil(size(SS_SavedData,1)*(trimData)):end,:);
            SS_accepted = SS_accepted(ceil(size(SS_accepted,1)*(trimData)):end,:);
        end
    
        f = figure; 
        hold on; grid off; xlabel("xd (m)"); ylabel("Force (N)"); xline(0, "--",'HandleVisibility','off'); yline(0, "--",'HandleVisibility','off');
        plot(SS_accepted(2:end,6), SS_SavedData(2:end,1), "DisplayName", "spring 1");  f.Position = [800 100 400 400];
        if dynamicSysParams.spring2_Active
            hold on; plot(-SS_accepted(2:end,6), SS_SavedData(2:end,2), "DisplayName", "spring 2");
        end
        legend show;
        legend('Location','northwest')
    end
else
    figure;
    t_ = ceil(size(odeParams.accepted_local,1)*.9);
    temp = [odeParams.accepted_local(t_:end,6), odeParams.accepted_local(t_:end,6)*dynamicSysParams.kd+ odeParams.accepted_local(t_:end,7)*dynamicSysParams.cd];
    plot(temp(:,1),temp(:,2))
    
end




