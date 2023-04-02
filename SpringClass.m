classdef SpringClass
    properties
        spring;
        name; % The name of the spring 
        d0; r0; % Wire diameter/Radius (mm)
        D0; R0; % Coil diameter/Radius (mm) 
        R; % The size of the coil in the midst of simulation
        l0; % Initial length of spring l0 (mm) 
        N; % Number of coils 
        alpha0; % Initial helix angle
        L; % The total length of the spring

        S_M_finish; % Martensite finish stress
        S_M_start; % Martensite start stress
        S_A_start; % Austenite finish stress
        S_A_finish; % Austenite start stress

        Gm; % Shear module of martensite
        Ga; % Shear module of austenite
        eL; % Max elastic strain of SMA material
        T; % Ambient temprature
        As; % Austenite start temprature
        Af; % Austenite finish temprature
        Ms; % Martensite start temprature
        Mf; % Martensite finish temprature
        
        z0; % Martensite ratio at start of the current phase
        zS0; % Martensite ratio induced by stress at start of the current phase (Equal to z if we only ahve superelastic movement)
        zT0; % Martensite ratio induced by temprature at start of the current phase (Equal to z if we only ahve superelastic movement)
        z; % Current martensite ratio of the spring crosssection
        elongation; % Current spring elongation
        strainDot; % Current force direction of the spring
        prevStrainDot; % Previous force direction of the spring
        phase; % Current phase of the spring
        previousePahse; % Previous phase of the spring
        n1; % Number of parts to divide the transformation part (Stress-wise)
        params; 
        verbose = false;
        Ka; % Stiffness of the spring when its full austenite
        Km; % Stiffness of the spring when its full martensite

        preStretchPercent; % The prestretch of the spring (In percentage of the range from 0 to S_M_Start)
        preStretchDisplacement; % The prestretch of in millimeters which will be calculated from preStretchPercent when we are initiating the spring
        stress; % The current stress in spring's cross section
        force;
        stressTable;
        strainTable;
        elong_SM_end; % The elongation at which the Austenite to Martensite transformation ends
        elong_SM_start; % The elongation at which the Austenite to Martensite transformation starts
        R_SM_end; % The coil diameter at spring_elong_SM_end

    end

    methods
        function obj = init(obj,d0,D0,l0,N,S_M_finish,S_M_start,S_A_start,S_A_finish,Gm,Ga,eL,T,As,Af,Ms,Mf,z0,zS0,zT0,z,eDot,phase,prevPhase,n1,preStretchPercent,params,name)
            if nargin ~= 0
                obj.name = name;
                obj.d0 = d0; obj.r0 = d0/2;
                obj.D0 = D0; obj.R0 = D0/2;
                obj.R = obj.R0;
                obj.l0 = l0;
                obj.N  = N ;
                obj.alpha0 = atan(obj.l0 / (2 * pi * obj.N * obj.R0) );
                obj.L = sqrt(obj.l0^2 + (2 * pi * obj.N * obj.R0)^2);
                obj.preStretchPercent = preStretchPercent;

                obj.S_M_finish = S_M_finish(1);
                obj.S_M_start = S_M_start(1);
                obj.S_A_start = S_A_start(1);
                obj.S_A_finish = S_A_finish(1);

                obj.Gm = Gm;
                obj.Ga = Ga;
                obj.eL = eL;
                obj.T = T;
                obj.As = As;
                obj.Af = Af;
                obj.Ms = Ms;
                obj.Mf = Mf;

                obj.z0 = z0; 
                obj.zS0 = zS0;
                obj.zT0 = zT0;
                obj.z = z; 
                obj.strainDot = eDot;
                obj.prevStrainDot = eDot;
                obj.phase = phase;
                obj.n1 = n1;
                obj.params = struct();
                obj.params.name = name;
                obj.params.Gm = Gm; 
                obj.params.Ga = Ga;
                obj.params.eL = eL;
                obj.params.z0 = z0;
                obj.params.zS0 = zS0;
                obj.params.S_M_finish = S_M_finish(1);
                obj.params.S_M_start = S_M_start(1);
                obj.params.S_A_finish = S_A_finish(1);
                obj.params.S_A_start = S_A_start(1);

                % Instantiating the constants necessary for finding the
                % stress when we are in transformation zone.
                q_A2M = linspace(S_M_finish(1),S_M_start(1),n1)';
                q_M2A = linspace(S_A_start(1),S_A_finish(1),n1)';
                [A2M_CORE, M2A_CORE] = martensiteRatio_CORE(obj, q_A2M, q_M2A);
                obj.stressTable = table(q_A2M,q_M2A,A2M_CORE,M2A_CORE);
                obj.strainTable = strainTableGenerator(obj, A2M_CORE, M2A_CORE, q_A2M, q_M2A, obj.params);

                % Calulating the elongation at which the Austenite
                % transformation ends and coil's radii when this happens.
                syms del; obj.elong_SM_end = vpasolve(obj.strainTable.strain_A2M(1) == del*obj.r0/2/pi/obj.N*3/4*cos(obj.alpha0)^2/obj.R0^2/cos(asin(del/obj.L+sin(obj.alpha0)))^2);
                obj.elong_SM_end = obj.elong_SM_end(obj.elong_SM_end>0);
                obj.elong_SM_end = double(obj.elong_SM_end(end));
                obj.R_SM_end = double(obj.R0 * cos(asin( obj.elong_SM_end / obj.L + sin(obj.alpha0) )) / cos(obj.alpha0));

                obj.elong_SM_start = vpasolve(obj.strainTable.strain_A2M(end) == del*obj.r0/2/pi/obj.N*3/4*cos(obj.alpha0)^2/obj.R0^2/cos(asin(del/obj.L+sin(obj.alpha0)))^2);
                obj.elong_SM_start = obj.elong_SM_start(obj.elong_SM_start>0);
                obj.elong_SM_start = double(obj.elong_SM_start(end));

                obj.preStretchDisplacement = preStretchPercent/100*obj.elong_SM_start * 0.001; % In meters

                obj.Ka = obj.d0^4 * obj.Ga / 8 / obj.D0^3 / obj.N * 1000;
                obj.Km = obj.d0^4 * obj.Gm / 8 / obj.D0^3 / obj.N * 1000;

                obj.params.z0 = z0;
                obj.params.zS0 = zS0;
                obj.params.forceLinear = params.forceLinear; % If we need to force the spring to act linear
                obj.params.forceLinearZ0 = params.forceLinearZ0; % The z0 just before we start forcing spring to be linear
                obj.params.forceLinearStress = params.forceLinearStress; % The stress just before we start forcing spring to be linear
                obj.params.forceLinearDir = params.forceLinearDir; % The direction at which the spring cant continiue linearly
                
                
                % Important: For the initial stress and force of the spring, we
                % approximate a 100% austenite spring with a prestretch
                obj.force = obj.preStretchDisplacement * obj.Ka;
                obj.stress = obj.force / 2 / pi / obj.r0^3 *3 *obj.R0;
            else
                disp("ERROR! No input for constructor")
            end
        end

        function [A2M_CORE, M2A_CORE] = martensiteRatio_CORE(obj,q_A2M,q_M2A)
            % This function only calculates the part of the formulas that doesn't
            % have z0 in them. Because the z0 changes on every phase changes, by
            % this we avoide redundant calculations. This function has to be
            % complemented by "strainTableGenerator" function.
            % ** Should run on application initialization.
        
            % Args:
            % q_A2M, q_M2A: matrix: The stresses at which austenite to martensite phase transform happens
        
            % Returns:
            % Two matrices containing the calculations
            
            % Defining the phase transformation as function handlers
            % 1.Martensite percentage for austenite to martensite
            Z_A2M_CORE = @(stress) cos(pi/(obj.S_M_start - obj.S_M_finish) * (stress - obj.S_M_finish));
            % 2.Martensite percentage for martensite to austenite
            Z_M2A_CORE = @(stress) (cos( pi/(obj.S_A_start-obj.S_A_finish)*(stress-obj.S_A_start) )+1);
            
            A2M_CORE = Z_A2M_CORE(q_A2M);
            M2A_CORE = Z_M2A_CORE(q_M2A);
        
        end

        function tbl =  strainTableGenerator(obj, A2M_CORE, M2A_CORE, q_A2M, q_M2A, params)
            % This function calculates the respective strain for each particular
            % stress.
            % ** Should run every time z0 chanegs.
        
            % Args:
            % A2M_CORE, M2A_CORE: matrix: Matrices that contain the core of the
            %   calculations. cos(pi/(S_M_start - S_M_finish) * (stress - S_M_finish))
            %   for austenite to martensite transform and cos(pi/(S_M_start - S_M_finish) * (stress - S_M_finish))
            %   for martensite to austenite.
            % q_A2M, q_M2A: double: The stresses at which austenite to martensite phase transform happens
        
            % Returns:
            % A table containig the respective strain for each stress in q_A2M or
            % q_M2A matrix.
            Z_A2M = (1-params.zS0)/2*A2M_CORE+(1+params.zS0)/2;
            strain_A2M = ((Z_A2M/params.Gm + (1-Z_A2M)/params.Ga)) .* q_A2M + params.eL * Z_A2M;
            
            Z_M2A = params.z0/2*M2A_CORE;
            strain_M2A = ((Z_M2A/params.Gm + (1-Z_M2A)/params.Ga)) .* q_M2A + params.eL * Z_M2A;
        
            tbl = table(strain_A2M, strain_M2A);
        end

    end

end