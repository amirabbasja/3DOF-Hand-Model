clear all -except adams;

system = sysParams();
system = system.init();

odeParams = struct(); % The variable that will be passed into ode15s
dynamicSysParams = struct(); % The variable that will be passed into dynamic model
stopParams = struct(); % The variable that will be passed into  ode15s but will be used in stopSim function

ampElbow = system.ampElbow;
freqElbow = system.freqElbow;
ampWrist = system.ampWrist;
freqWrist= system.freqWrist;

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
dynamicSysParams.freqElbow = system.freqElbow;
dynamicSysParams.freqWrist = system.freqWrist;

% Necessary for not running into errors
dynamicSysParams.spring1_Active = false;
dynamicSysParams.spring2_Active = false;

%%%%%% Override data
dynamicSysParams.ampWrist = 1;
dynamicSysParams.ampElbow = 0;
dynamicSysParams.freqElbow = 5;
dynamicSysParams.freqWrist = 5;
%%%%%%

% The initial state of the system
th1_Init = 0;
Dth1_Init = 0;
th2_Init = 0;
Dth2_Init = 0;

% Defining the dummy spring global variables. we instantiate these
% variables so that we can use ode15s (modified version) 
global hasSMASpring;
hasSMASpring = false;

% Defining the matrices that will save the iteration datas
global accepted savedData;
savedData = zeros(1,1);
accepted = zeros(1,5);

% Defining the parameters necessary form event fucntion to stop the
% integration
peakRepeat_th1 = 4;
peakRepeat_th2 = 4;
stopParams.peakDelay = 1; % After the stopping criterion is met, we delay the stopping of the system to be more sure of the steadiness of the amplitude. 
stopParams.peakTollerance = 0.005; % The amount of peak reduction relative to its previous peak to stop the symulation
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

tStart = 0;
tEnd = 320;
increment = 0.001;

odeOptions = odeset('RelTol', 1e-3, 'AbsTol', 1e-5);

% Locating the solver folder and adding it to path
addpath(sprintf("./solver/R%s/", "2019b"))
tic
[odeParams,t,Y] = ode15s_Modified( ...
    @(t,Y,params) dynamicEqn_2DOF_DEAD_MASS(t,Y,params), ...
    linspace(tStart,tEnd,(tEnd-tStart)/increment), ...
    [th2_Init;Dth2_Init;th1_Init;Dth1_Init], odeOptions, dynamicSysParams, odeParams, dynamicSysParams);
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
disp(steadyStateAmp_th1)
out = [];
out.tout = t; 
out.th2 = Y(:,1); 
out.Dth2= Y(:,2); 
out.th1 = Y(:,3); 
out.Dth1 = Y(:,4); 

hold on; grid off; plot(out.tout, out.th2); 
