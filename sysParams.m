% The design parameters of the system that can be fed into the SIMULINK
% simulation or matlab scripts

classdef sysParams
    properties
        %%% Hand
        m1; 
        m2; 
        g; 
        c1; 
        c2; 
        kt1Mul; 
        kt2Mul; 
        l1; 
        l2; 
        lc1; 
        lc2; 
        ld;
        ld2;
        xdInit; 
        i1; 
        i2; 

        cd; 
        kd; 
        md; 
        
        cd2;
        kd2;
        md2;

        ampElbow;
        ampWrist;
        freqElbow;
        freqWrist;

        freqSweepe;

        % SMA Springs: Each sprig is a struct, consisted of the following
        % properties:
        % d0; % r0; % Wire diameter/Radius (mm)
        % D0; % R0; % Coil diameter/Radius (mm) 
        % l0; % Initial length of spring l0 (mm) 
        % N; % Number of coils 
        % alpha0; % Initial helix angle
        % L; % The total length of the spring
        % Mf; % Ms; % As; % Af;
        % CA;
        % CM ;
        % critStressFinish;
        % critStressStart;
        % Ga;
        % Gm;
        % eL;
        % T;
        % S_M_finish;
        % S_M_start;
        % S_A_start;
        % S_A_finish; 
        % preStretchDisplacement;
        % preStretchPercentage;
        % Ka;
        % Km;
        % Active; % If the system uses spring1, true, else false.
        % elong_SM_end; % The elongation at which the Austenite to Martensite transformation changes
        % R_SM_end; % The coil diameter at spring_elong_SM_end

        % Spring1
        spring1;
        % Spring2
        spring2;
    end
    
    methods
        function obj = init(obj)
            obj.m1 =1.1826; 
            obj.m2 =0.446; 
            obj.g =9.80665; 
            obj.c1 =0.1; 
            obj.c2 =0.02; 
            obj.kt1Mul =2; 
            obj.kt2Mul =2; 
            obj.l1 =0.2513; 
            obj.l2 =0.1899; 
            obj.lc1 =0.1166; 
            obj.lc2 =0.0701; 
            obj.ld =.1; 
            obj.ld2 = .1;
            obj.xdInit =0.0; 
            obj.i1 =5824.9e-6; 
            obj.i2 =1304.9e-6; 

            obj.cd =0; 
            obj.kd =0; 
            obj.md =0;

            obj.kd2 = 50;
            obj.cd2 = 1;
            obj.md2 = 0.15;
            
            obj.ampElbow = 0;
            obj.freqElbow = 0;
            obj.ampWrist = 0;
            obj.freqWrist = 0;

            obj.freqSweepe = [0 1 1.5 2 2.5 3 3.5 4 4.2 4.4 4.5 4.6 4.8 4.9 5 5.1 5.3 5.4 5.5 5.6 5.75 6 6.1 6.2 6.3 6.4 6.5 6.6 6.75 6.9 7 7.1 7.2 7.25 7.35 7.5 7.6 7.8 8 8.1 8.3 8.5 8.6 8.75 8.9 9 9.1 9.2 9.25 9.3 9.4 9.5 9.6 9.7 9.8 9.9 10 10.5 11 11.5 12 12.5 13 13.5 14 14.5 15 15.5 16 ];


            %% Spring 1
            obj.spring1 = struct();
            obj.spring1.name = "spring1";
            obj.spring1.d0 = ;  obj.spring1.r0 = obj.spring1.d0/2; % Wire diameter/Radius (mm)
            obj.spring1.D0 = ;  obj.spring1.R0 = obj.spring1.D0/2; % Coil diameter/Radius (mm) 
            obj.spring1.l0 = ; % Initial length of spring l0 (mm) 
            obj.spring1.N = ; % Number of coils 
            obj.spring1.alpha0 = atan(obj.spring1.l0 / (2 * pi * obj.spring1.N * obj.spring1.R0) ); % Initial helix angle
            obj.spring1.L = sqrt(obj.spring1.l0^2 + (2 * pi * obj.spring1.N * obj.spring1.R0)^2); % The total length of the spring

            % No Mf is provided butit has very little effect on the diagrams so we
            % extrapolate and assign 0 to it
            obj.spring1.Mf = ; obj.spring1.Ms = ; obj.spring1.As = ; obj.spring1.Af = ;
            
            obj.spring1.CA = ;
            obj.spring1.CM = ;
            obj.spring1.critStressFinish = ; obj.spring1.critStressStart =  ;
            
            obj.spring1.Ga = ; obj.spring1.Gm = ; obj.spring1.eL = ;
            obj.spring1.T = ;
            
            % Coefficients to fit the analytical data to experimental data
            % Done by myself, Not the paper
            b = [1 1 1 1];

            obj.spring1.S_M_finish = obj.spring1.critStressFinish + obj.spring1.CM*(obj.spring1.T-obj.spring1.Ms) * b(1);
            obj.spring1.S_M_start = obj.spring1.critStressStart + obj.spring1.CM*(obj.spring1.T-obj.spring1.Ms) * b(2);
            obj.spring1.S_A_start = obj.spring1.CA*(obj.spring1.T-obj.spring1.As) * b(3);
            obj.spring1.S_A_finish = obj.spring1.CA*(obj.spring1.T-obj.spring1.Af) * b(4); 

            obj.spring1.Ka = obj.spring1.d0^4 * obj.spring1.Ga / 8 / obj.spring1.D0^3 / obj.spring1.N * 1000; % Stiffness of the spring when its full austenite
            obj.spring1.Km = obj.spring1.d0^4 * obj.spring1.Gm / 8 / obj.spring1.D0^3 / obj.spring1.N * 1000; % Stiffness of the spring when its full martensite

            obj.spring1.preStretchPercentage = ; % In percentage; 0 for 0Mpa and 100% for S_M_start (Mpa). As you notice, we have taken the initial state of the spring to be full austenite
            obj.spring1.Active = false; % Default is that the system doesn't have SMA spring
            
            if obj.spring1.Active; fprintf("Spring 1 approx. Stiffness: (A) %.2f || (M) %.2f [N/m]\n", obj.spring1.Ka,  obj.spring1.Km); end
            
            %% Spring 2
            % We have taken spring1 and spring2 to be identical
            obj.spring2 = obj.spring1;
            obj.spring2.name = "spring2";

            if obj.spring2.Active; fprintf("Spring 2 approx. Stiffness: (A) %.2f || (M) %.2f [N/m]\n", obj.spring2.Ka,  obj.spring2.Km); end
            %% Override data
        end
    end
end