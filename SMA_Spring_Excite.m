function [excitedSpring,stress,force] = SMA_Spring_Excite(springObj, Elongation)
    % Only account for change of coil diameter if the stress is below
    % the transfromations (A->M) stress. This is to avoid exponential
    % increase of stress after transformation.
    if abs(Elongation) <= abs(springObj.elong_SM_end)
        % Calculating the spring diameter change
        alpha = asin( abs(Elongation) / springObj.L + sin(springObj.alpha0) );
        R = springObj.R0 * cos(alpha) / cos(springObj.alpha0);
    else
        R = springObj.R_SM_end;
    end

    % The effective strain is equal to (3/4) of maximum shear strain
    shearStrain = Elongation * springObj.r0 / (2 * pi * springObj.N * R^2) * (3/4);
   
    if abs(springObj.elongation) <= abs(Elongation)
        eDot = 1;
        springObj.strainDot = eDot;
    else
        eDot = -1;
        springObj.strainDot = eDot;
    end

    [z0,z,phase,stress,strainTable,params] = MAIN( ...
        springObj.z0, ...
        springObj.z, ...
        springObj.phase, ...
        shearStrain, ...
        eDot, ...
        springObj.prevStrainDot, ...
        springObj.stress, ...
        springObj.T, springObj.eL, springObj.Gm, springObj.Ga, springObj.Ms, springObj.As, springObj.S_M_finish, springObj.S_M_start, springObj.S_A_finish, springObj.S_A_start, springObj.verbose, springObj.stressTable, springObj.strainTable, springObj.params ...
        ,Elongation);

    % Note that because the force that the spring excerts is in the
    % opposite direction of excitation, we have to multiply the result
    % by -1 
    force = -stress * 2 * pi * springObj.r0 ^ 3 / 3 / R;
    
    excitedSpring = springObj;
    excitedSpring.stress = stress;
    excitedSpring.force = force;
    excitedSpring.z0 = z0;
    excitedSpring.z = z;
    excitedSpring.phase = phase;
    excitedSpring.strainTable = strainTable;
    excitedSpring.prevStrainDot = eDot;
    excitedSpring.elongation = Elongation;
    excitedSpring.R = R;
    excitedSpring.params = params;

    if springObj.verbose
            fprintf("%s : Elongation: %.5f  |  stress: %.2f  |  strainDot: %.f | phase: %.f  |  z: %.5f  | z0: %.5f | shearStrain: %.8f\n_______________________________________________________________\n", ...
                springObj.name ,Elongation, stress, springObj.strainDot, phase, z, z0, shearStrain);
    end

    function [z0,z,phase,stress,strainTable,params] = MAIN(z0,z,previousePahse,strain,strainDot, prevStrainDot, prevStress,T,eL,Gm,Ga,Ms,As, S_M_finish, S_M_start, S_A_finish, S_A_start, verbose, stressTable, strainTable, params, elongation)
        % This function with strain as an input and some other inputs, returns
        % stress and some other variables (These variables are used as inputs for the next iteration)
        % This function helps us avoide doing the calculations inside the loop
        % INPUTS:
        % z0: The martensite ratio at start of previouse phase
        % z: The martensite ratio of previouse iteration (It is different form z0)
        % prevPhase: The phase of previouse iteration
        % strain: Current strain
        % strainDot: Current strain rate
    
        % Only pseudoelastic
        zT = 0; zT0 = 0;
        params.z = z;
    
        % When we reverse in the middle of phase transformation, we need to
        % update z0, but due to us not having the stress, we cant know if we
        % have changed phase (phaseFinder function needs stress as an input)
        % to solve this error, we know if we have reversed if previose eDot
        % is not equal to current eDot. when they are not equal, we update
        % the value of z0
        reverse = false;
        if prevStrainDot * strainDot < 0
            reverse = true;
        end

        if reverse
            if previousePahse == 2 || previousePahse == 5
                params.forceLinearZ0 = z0;
                params.forceLinearStress = abs(prevStress);
            end
            
            if previousePahse ~= 2 && previousePahse ~= 5 && previousePahse ~= 1 && previousePahse ~= 6
                % Note that when z == 1 or 0, and we have reversion, we don't
                % need to force the spring to act linear, it already is.
                if 0 < strainDot && S_A_start(1) < abs(prevStress) && params.forceLinear == false && (z ~= 0 && z ~= 1)
                    params.forceLinear = true;
                    params.forceLinearDir = 1;
                elseif strainDot < 0 && abs(prevStress) < S_M_start(1) && params.forceLinear == false && (z ~= 0 && z ~= 1)
                    params.forceLinear = true;
                    params.forceLinearDir = -1;
                else
                    params.forceLinear = false;
                    params.forceLinearDir = nan;
                end
            end
    
            zS0 = z;
            z0 = z;
            params.z0 = z0;
            params.zS0 = zS0;
            strainTable = strainTableGenerator(stressTable.A2M_CORE, stressTable.M2A_CORE, stressTable.q_A2M, stressTable.q_M2A, params);
        end

        [stress, z0, params] = stressFinder(previousePahse, strain, strainDot, prevStress, T, eL, Gm, Ga, z0, zT, S_M_finish, S_M_start, S_A_finish, S_A_start, verbose, stressTable, strainTable, params, elongation);

        phase = phaseFinder(previousePahse,  stress, strainDot,z,z0, S_M_finish, S_M_start, S_A_finish, S_A_start, verbose, params);
    
        if previousePahse ~= phase
            phaseChange = true;
        else
            phaseChange = false;
        end
        
        % When we finsh phase 2 and go to phase 3 z is not exactly equal to
        % 1, this causes errors in stressFinder function, so by this line we
        % fix this
        if ((previousePahse == 2 || previousePahse == 21) && phase == 3)
            z = 1;
        elseif ((previousePahse == 5 || previousePahse == 51) && phase == 6) 
            z = 0;
        end

        if phase == 8000
            if params.forceLinearDir == strainDot
                z0 = params.forceLinearZ0;
                zS0 = params.forceLinearZ0;
                params.z0 = params.forceLinearZ0;
                params.zS0 = params.forceLinearZ0;
                strainTable = strainTableGenerator(stressTable.A2M_CORE, stressTable.M2A_CORE, stressTable.q_A2M, stressTable.q_M2A, params);
            end
            params.forceLinear = false;
            params.forceLinearStress = nan;
            params.forceLinearZ0 = nan;
            params.forceLinearDir = nan;
            
        end

        if phaseChange && phase ~= 8000 && previousePahse ~= 8000 && phase ~= 800 && previousePahse ~= 800
            zS0 = z;
            z0 = z;
            params.z0 = z0;
            params.zS0 = zS0;
            strainTable = strainTableGenerator(stressTable.A2M_CORE, stressTable.M2A_CORE, stressTable.q_A2M, stressTable.q_M2A, params);
        end
    
        %The get_z function cant process negative stresses, so we multiply the
        %stress by -1 and at the end we change it back
        if(strain < 0)
            NEGATIVESTRAIN = true;
            stress = stress * -1;
        end

        z = get_z(stress,strainDot,T,Ms,As,z0,zT,z, S_M_finish, S_M_start, S_A_finish, S_A_start, params);

        if(strain < 0)
            NEGATIVESTRAIN = true;
            stress = stress * -1;
        end
    end

    function [stress, z0, params] = stressFinder(prevPhase, strain, strainDot, prevStress, T, eL, Gm, Ga, z0, zT0, S_M_finish, S_M_start, S_A_finish, S_A_start, verbose, stressTable, strainTable, params, elongation)
        % In this Function we try to get the stress with strain as an input. the
        % logic is: First we assume we are not in transformation phase and try
        % to solve the equation, if no answers were found, our assumption was
        % wrong and we are in transformation phase
        
        % z0 is the z at the start of new phase, but z is the previouse step's
        % martnsite ratio. this comes handy when we are decreasing the strain
        zS0 = z0 - zT0;
    
        % The logic for solving the equations: When we are in the linear parts
        % of the loading process, we avoide using vpasolve and in the new
        % update, we substitute the parameters in the stress equation and solve
        % for stress. For the non-linear parts, Again, we have avodied using
        % vpasolve. We divide the transformation regions(A->M and m->A) in to 
        % equal parts and then we finde the respective strain for each stress.
        % Then when we want to solve the transformation stress for each strain, 
        % we use the data acquired previousely. This causes the app's speed to
        % increas like 20 times without any loss of generality.
    
        % UPDATE: some times, we have negative strain. In order for this
        % function to handle that, we define This constant. if it is true, all
        % the calculations will be in positive numbers but at the end they will
        % be multiplied by -1
        NEGATIVESTRAIN = false;
        stress = nan; % Pre assigning the stress
        
        if(strain < 0)
            NEGATIVESTRAIN = true;
            strain = strain * -1;
        end

        if params.forceLinear
            s = Func_SMA_Shear_Modulus(z0,Gm,Ga,1) * (strain - eL*z0);
            if params.forceLinearStress < s && S_M_start < prevStress && 0 < strainDot && size(s,1) == 0
                % The linear stress-strain line shouldn't cross the
                % transformation, if we reach there, we should use 
                % transformatio formulas
                temp_params = params;
                temp_params.zS0 = params.forceLinearZ0;
                temp_params.z0 = params.forceLinearZ0;
                tempStrainTable = strainTableGenerator(stressTable.A2M_CORE, stressTable.M2A_CORE, stressTable.q_A2M, stressTable.q_M2A, temp_params);
                s = stressFromStrainArray(strain, strainDot, tempStrainTable.strain_A2M, tempStrainTable.strain_M2A, stressTable.q_A2M, stressTable.q_M2A);
    
                % Some times z0 is very small or very close to 1, and we are
                % coming from phase 51 to 6 or 21 to 3. this causes 
                % application to be lost because it didnt see phase
                % 8000 to stop forcing linearization of spring. we solve this here
                if size(s,1) == 0 
                    r_ = Func_SMA_Shear_Modulus(1,Gm,Ga,1) * (strain - eL*z0);
                    if (S_M_finish <= r_) && (r_ <= inf); s = [r_]; else s = []; end
                end
            elseif s < params.forceLinearStress && prevStress < S_A_start && strainDot < 0 && size(s,1) == 0
                % The linear stress-strain line shouldn't cross the
                % transformation, if we reach there, we should use
                % transformatio formulas
                temp_params = params;
                temp_params.zS0 = params.forceLinearZ0;
                temp_params.z0 = params.forceLinearZ0;
                tempStrainTable = strainTableGenerator(stressTable.A2M_CORE, stressTable.M2A_CORE, stressTable.q_A2M, stressTable.q_M2A, temp_params);
                s = stressFromStrainArray(strain, strainDot, tempStrainTable.strain_A2M, tempStrainTable.strain_M2A, stressTable.q_A2M, stressTable.q_M2A);
    
                % Some times z0 is very small or very close to 1, and we are
                % coming from phase 51 to 6 or 21 to 3. this causes 
                % application to be lost because it didnt see phase
                % 8000 to stop forcing linearization of spring. we solve this here
                if size(s,1) == 0 
                    r_ = Func_SMA_Shear_Modulus(0,Gm,Ga,1) * (strain - eL*z0);
                    if (0 <= r_) && (r_ <= S_A_finish); s = [r_]; else s = []; end
                end
            end
    
            if size(s,1) == 1
                stress = s;
            end
        end
    
        if 0 <= strainDot && ~params.forceLinear
            % If stress/strain is increasing search like this
            % 1. Assume no transformation
            if z0 == 0
                r_ = Func_SMA_Shear_Modulus(z0,Gm,Ga,1) * (strain - eL*z0);
                if (0 <= r_) && (r_ <= S_M_start); s = [r_]; else s = []; end
            elseif z0 == 1
                r_ = Func_SMA_Shear_Modulus(z0,Gm,Ga,1) * (strain - eL*z0);
                if (S_A_start <= r_) && (r_ <= inf); s = [r_]; else s = []; end
            elseif 0 < z0 && z0 < 1
                r_ = Func_SMA_Shear_Modulus(z0,Gm,Ga,1) * (strain - eL*z0);
                if (0 <= r_) && (r_ <= S_M_start); s = [r_]; else s = []; end
            end
    
            if(size(s,1) == 0)
                %2.No answers found, Try phase transformation
        
                % **Note**: I dont know why, but when we use (Func_SMA_Shear_Modulus(Z_A2M(q),Gm,Ga,1)^-1)
                % instead of (Z_M2A(q)/Gm + (1-Z_M2A(q))/Ga) in the relation below, vpasolve has some
                % convergance issues. 
                s = stressFromStrainArray(strain, strainDot, strainTable.strain_A2M, strainTable.strain_M2A, stressTable.q_A2M, stressTable.q_M2A);
    
                if(size(s,1) == 0)
                    % When we are at the very first step of the transformation,
                    % we dont know if the z0 has changed or not, so first we
                    % try with the z0 that we had (previouse if statement) and
                    % if we get no answers from it, we try assigning the
                    % current z to z0 and try the transformation equation
                    % again.
                    temp_params = params;
                    temp_params.zS0 = params.forceLinearZ0;
                    temp_params.z0 = params.forceLinearZ0;
                    tempStrainTable = strainTableGenerator(stressTable.A2M_CORE, stressTable.M2A_CORE, stressTable.q_A2M, stressTable.q_M2A, temp_params);
                    s = stressFromStrainArray(strain, strainDot, tempStrainTable.strain_A2M, tempStrainTable.strain_M2A, stressTable.q_A2M, stressTable.q_M2A);
                    if (size(s,1) == 0)
                        if(prevPhase == 2 || prevPhase == 21 || prevPhase == 8000)
                            % No answers found. Some times this occurs when we finish phase
                            % transformation and continue lineraly in the same
                            % direction. At the point of transition two problems occur. 
                            % first z never reaches exactly to 1 (despite how close it may be) 
                            % ,so we cant really know when the transformation ends. 
                            % but when iteration  through e(i)'s we reach a point where no
                            % answer is found from vpasolve. This is where we know we
                            % have passed the transformation level. second is that we
                            % have to update z0 to continue iteration, but we have to
                            % do it without breaking the iteration loop. so we use
                            % recursion to solve this matter.
                            stress = stressFinder(prevPhase, strain, strainDot, prevStress, T, eL, Gm, Ga, 1, zT0, S_M_finish, S_M_start, S_A_finish, S_A_start, verbose, stressTable, strainTable,params);

                            if size(stress,1) == 1
                                z0 = 1;
                                params.z0 = 1;
                                params.zS0 = 1;
                            end
                        else
                            if verbose; disp("No answers found! Transformation"); end
                            if verbose; fprintf("strain = %.5f  || z0 = %.5f || prevPhase = %.1f|| strainDot = %.1f\n",strain,z0,prevPhase,strainDot); end
                        end
                    else
                        stress = s(1);
                        if verbose; disp("One answer found! The very first step of Transformation"); end
                    end
            
                elseif(size(s,1) == 1)
                    %Only one answer found
                    stress = s(1);
                    if verbose;  disp("One answer found! Transformation"); end
                else
                    %More than one answer found! 
                    %Print an Error message!
                    if verbose; disp("#########ERROR##########\nMore than one answer found! Transformation"); end
                    stress = -1;
                end
        
            elseif(size(s,1) == 1)
                %Only one answer found
                stress = s(1);
            else
                %More than one answer found! 
                %Print an Error message!
                if verbose; disp("#########ERROR##########\nMore than one answer found!"); end
                stress = -1;
            end
    
    
        
        elseif strainDot < 0 && ~params.forceLinear
            % If stress/strain is decreasing search like this
            % 1.Assume no transformation
            if z0 == 0

                r_ = Func_SMA_Shear_Modulus(z0,Gm,Ga,1) * (strain - eL*z0);
                if (0 <= r_) && (r_ <= S_A_finish); s = [r_]; else s = []; end
    
                % When the wire turns back at phase 1 with z=0; the code would
                % think its in transformation phase. we add this line to solve
                % this problem. It has been added as a quick fix and its stability 
                % has to be checked more. 
                if prevPhase == 1 || prevPhase == 11
                    r_ = Func_SMA_Shear_Modulus(z0,Gm,Ga,1) * (strain - eL*z0);
                    if (0 <= r_) && (r_ <= S_M_start); s = [r_]; else s = []; end
                end
            elseif z0 == 1
                r_ = Func_SMA_Shear_Modulus(z0,Gm,Ga,1) * (strain - eL*z0);
                if (S_A_start <= r_) && (r_ <= inf); s = [r_]; else s = []; end
            elseif 0 < z0 && z0 < 1
                r_ = Func_SMA_Shear_Modulus(z0,Gm,Ga,1) * (strain - eL*z0);
                if (S_A_start <= r_) && (r_ <= inf); s = [r_]; else s = []; end
            end
        
            if(size(s,1) == 0)
                % 2.No answers found, Try phase transformation
        
                % **Note**: I dont know why, but when we use (Func_SMA_Shear_Modulus(Z_A2M(q),Gm,Ga,1)^-1)
                % instead of (Z_M2A(q)/Gm + (1-Z_M2A(q))/Ga) in the relation below, vpasolve has some
                % convergance issues. 
                s = stressFromStrainArray(strain, strainDot, strainTable.strain_A2M, strainTable.strain_M2A, stressTable.q_A2M, stressTable.q_M2A);
    
                if(size(s,1) == 0)
                    % When we are at the very first step of the transformation,
                    % we dont know if the z0 has changed or not, so first we
                    % try with the z0 that we had (previouse if statement) and
                    % if we get no answers from it, we try assigning the
                    % current z to z0 and try the transformation equation
                    % again.
                    temp_params = params;
                    temp_params.zS0 = params.forceLinearZ0;
                    temp_params.z0 = params.forceLinearZ0;
                    tempStrainTable = strainTableGenerator(stressTable.A2M_CORE, stressTable.M2A_CORE, stressTable.q_A2M, stressTable.q_M2A, temp_params);
                    s = stressFromStrainArray(strain, strainDot, tempStrainTable.strain_A2M, tempStrainTable.strain_M2A, stressTable.q_A2M, stressTable.q_M2A);
    
                    if (size(s,1) == 0)
                        if(prevPhase == 5 || prevPhase == 51 || prevPhase == 8000)
                            % No answers found. Some times this occurs when we finish phase
                            %  transformation and continue lineraly in the same
                            % direction. At the point of transition two problems occur. 
                            %  first z never reaches exactly to 1 (despite how close it may be) 
                            % ,so we cant really know when the transformation ends. 
                            %  but when iteration  through e(i)'s we reach a point where no
                            % answer is found from vpasolve. This is where we know we
                            % have passed the transformation level. second is that we
                            % have to update z0 to continue iteration, but we have to
                            % do it without breaking the iteration loop. so we use
                            % recursion to solve this matter.
                            stress = stressFinder(prevPhase, strain, strainDot, prevStress, T, eL, Gm, Ga, 0, zT0, S_M_finish, S_M_start, S_A_finish, S_A_start, verbose, stressTable, strainTable,params);
                            
                            if size(stress,1) == 1
                                z0 = 0;
                                params.z0 = 0;
                                params.zS0 = 0;
                            end
                        else
                            if verbose; disp("No answers found! Transformation"); end
                            if verbose; fprintf("strain = %.5f  || z0 = %.5f || prevPhase = %.1f|| strainDot = %.1f\n",strain,z0,prevPhase,strainDot); end
                        end
                    else
                        stress = s(1);
                        if verbose; disp("One answer found! The very first step of Transformation"); end
                    end
            
                elseif(size(s,1) == 1)
                    %Only one answer found
                    stress = s(1);
                    if verbose; disp("One answer found! Transformation"); end
                else
                    %More than one answer found! 
                    %Print an Error message!
                    if verbose; disp("#########ERROR##########\nMore than one answer found! Transformation"); end
                    stress = -1;
                end
        
            elseif(size(s,1) == 1)
                %Only one answer found
                stress = s(1);
            else
                %More than one answer found! 
                %Print an Error message!
                if verbose; disp("#########ERROR##########\nMore than one answer found!"); end
                stress = -1;
            end
        end
        
        if NEGATIVESTRAIN && ~isnan(stress)
            stress = stress * -1;
        end

        if isnan(stress)
            error("No stress was found in %s The parameters of the simulation are as follows:\n" + ...
                "Elongation %.5f\nprevPhase: %d \nstrain %.7f \nstrainDot %d \nprevStress %.3f\n\nparams: \n" + ...
                "z0   %5f\nforceLinear   %d\nforceLinearStress   %.5f\nforceLinearDir   %d\nforceLinearZ0   %.5f\n", params.name, elongation, prevPhase, strain, strainDot, prevStress,params.z0, params.forceLinear, params.forceLinearStress, params.forceLinearDir, params.forceLinearZ0)
        end
    
    end
    
    function [phase] = phaseFinder(prevPhase, stress, strainDot,z,z0, S_M_finish, S_M_start, S_A_finish, S_A_start, verbose, params)
        % Note that due to nature of my project zT is always equal to zero and we can use z instead of zS
    
        % This function returns the SMA's phase number in pseudoelastic a
        % procedure with Brinsons's (1993 and 1996) constitutive modeling
        % Phase Numbering:
        % loading before phase change : 1  || loading with phase change : 2  || loading after phase change : 3 
        % unloading before phase change : 4  || unloading with phase change : 5  || unloading after phase change : 6 
        % Phases for handeling inner loops are:
        % Loading while in middle of Martensite->Austenite phase change: 8
        % Loading when below martensit transformation start stress while  z != 0: 9
        % Loading while coming back into austenite phase Chaneg area while  z !=0: 10
        % **Note that the phases that are caused from reversals in loading in the middle of  a phase are repeated numers (e.g. 77 or 11)
        % Unloading before Austenite has changed to martensite (reversing in middle of phase 1):11
        % Unloading while in middle of Austenite->Martensite phase change: 22
        % Reversing phase 7 and continuing loading:77
        % Reversing in middle of phase 5 and continuing loading:55
    
    
        % For defining the inner loops we also need to know where we came from,
        % so we need the previous step's direction (loading or unloading)
        % * Note that stressDot is the direction of loading, +1 for loading and
        % -1 for unloading
    
        % UPDATE: some times, we have negative strain. In order for this
        % function to handle that, we define This constant. if it is true, all
        % the calculations will be in positive numbers but at the end they will
        % be multiplied by -1
        NEGATIVESTRAIN = false;
        
        if(stress < 0)
            NEGATIVESTRAIN = true;
            stress = stress * -1;
        end
    
        phase = -1;
    %     if (S_M_start < stress && stress < S_M_finish) && 0 < strainDot && z0 == 0 && params.forceLinear == false
    %         phase = 2;
        if (S_M_start < stress && stress < S_M_finish) && 0 < strainDot && params.forceLinear == false
            phase = 2;
        elseif not(S_M_start < stress && stress < S_M_finish) && (stress < S_M_start) && 0 < strainDot && z == 0
            phase = 1;
        elseif not(S_M_start < stress && stress < S_M_finish) && (S_M_finish < stress) && 0 < strainDot
            phase = 3;
        elseif not( S_A_finish < stress && stress < S_A_start  ) && strainDot < 0 && (S_A_start < stress)  && z == 1
            phase = 4;
        elseif ( S_A_finish < stress && stress < S_A_start  ) && strainDot < 0 && 0 < z && params.forceLinear == false
            phase = 5;
        elseif not( S_A_finish < stress && stress < S_A_start  ) && strainDot < 0  && (stress < S_A_finish) && prevPhase ~= 11
            %Some times we start unloading while in phase 1, so to
            %differenciate this with phase 6, we have to make a new phase and
            %"&& prevPhase ~= 1" condition as well. to avoid devision by zero
            %add
            phase = 6;
        elseif not(S_M_start < stress && stress < S_M_finish) && (stress < S_M_start) && strainDot < 0 && z == 0
            %Unloading before Austenite has changed to martensite (reversing in middle of phase 1)
            phase = 11;
        elseif (z0 < 1) && (stress < S_M_finish) && (S_A_start < stress) && strainDot < 0
            %Inner loop: Unloading while in middle of Austenite->Martensite phase change
            phase = 22;
        elseif ( S_A_finish < stress && stress < S_A_start  ) && 0 < strainDot && 0 < z
            phase = 55;
    %     elseif (0 < z0) && (z0 < 1) && (stress < S_M_finish) && (S_A_start < stress) && 0 < strainDot && params.forceLinear == false
    %         % When we are in the transformation region with z0 ~= 0 and we are
    %         % not forcing the spring to act linear
    %         phase = 222;
        elseif (0 < z0) && (z0 < 1) && (stress < S_M_finish) && (S_A_start < stress) && 0 < strainDot && params.forceLinear == true && stress <= params.forceLinearStress
            % When we are in the transformation region with z0 ~= 0 and we are
            % forcing the spring to act linear
            phase = 21;
        elseif (0 < z0) && (z0 < 1) && (stress < S_M_start) && (S_A_finish < stress) && strainDot < 0 && params.forceLinear == true && params.forceLinearStress <= stress
            % When we are in the transformation region with z0 ~= 0 and we are
            % forcing the spring to act linear
            phase = 51;
        elseif ((params.forceLinearStress < stress && 0 < strainDot && S_M_start < stress) ...
                || (stress < params.forceLinearStress && strainDot < 0 && stress < S_A_start && strainDot == params.forceLinearDir)) ...
                && params.forceLinear == true
            % When we have forced spring to act linear and we have passed the
            % previouse transformation stress. We have to change the parameters
            % of the spring to not force linear. 
            phase = 8000;
        elseif ( S_A_finish < stress && stress < S_A_start  ) && strainDot < 0 && z < 1
            %Inner loop: Unloading Phase transform from a point where z < 1
            phase = 8;
        elseif (S_M_start < stress && stress < S_M_finish) && 0 < strainDot && 0 < z
            %Inner loop: Loading while coming back into austenite phase Chaneg area while  z !=0
            phase = 10;
        elseif (S_A_start < stress && stress < S_M_start) && S_A_start < S_M_start
            % The stress region between S_M_start and S_A_start which no
            % transformation happens, Note that S_A_start has to be smaller
            % than S_M_start
            phase = 800;
        end
    
        if verbose; fprintf("Phase: %d\n", phase); end
    end

    function [z,zS,zT] = Func_Brinson_Conversion_To_Austenite(T,stress,As,zS0,zT0,currZ, S_M_finish, S_M_start, S_A_finish, S_A_start, params)
        % This function returns the amount of stress induced and tempreature induced
        % Martensite in the process of SMA's conversion to austenite.
        % The brinson's 1993 constitutive relation has been used
        % z0 : martensite ration in the beginning of loading cycle
        % zT0 : temprature induced martensite ration in the beginning of loading cycle
        % zS0 : Stress induced martensite ration in the beginning of loading cycle
        % currZ: Current martensite ratio of the sample; thiw helps us to avoid
        % division by zero in phase 11
        % Note that Af < T
        % Note that this relation stands true only if: CA*(T-As) < stress < CA*(T-Af)
        % Also note that phase transformation from Martensite to Austenite
        % happens only when specimen is in the "Unloading" phase and stress is
        % decreasing (zDot < 0)
    
        % Update: in this update we removed the criticalStressStart and
        % criticalStressFinish from the required arguments. the user has to
        % only provide S_A_start and S_A_finish (Critical stress that
        % conversion to austenite starts and ends respectivelt). the user has
        % to calculate these parameters and pass them to method. In some papers
        % there is no criticalStress provided (unlike the brinson's model. we use
        % this method to make the function more friendly for implementing other
        % papaers)
    
        % Update: The function will nolonger update z when we are forcing the
        % spring to act linear
        if (As < T) && (S_A_finish < stress && stress < S_A_start ) && ~params.forceLinear
            if 0 < currZ
                z0 = zT0 + zS0;
                z = z0/2*(cos( pi/(S_A_start-S_A_finish)*(stress-S_A_start) )+1);
                zS = zS0 - zS0/z0*(z0-z);
                zT = zT0 - zT0/z0*(z0-z);
            else
                zS = zS0;
                zT = zT0;
                z = zS + zT;
            end
            
        else
            zS = zS0;
            zT = zT0;
            z = zS + zT;
        end
    end
    
    function [z,zS,zT] = Func_Brinson_Conversion_To_Detwinned_Martensite(T,stress,Ms,zS0,zT0, S_M_finish, S_M_start, S_A_finish, S_A_start, params)
        % This function returns the amount of stress induced and tempreature induced
        % Martensite in the process of SMA's conversion to detwinned martensite.
        % The brinson's 1993 constitutive relation has been used
        % z0 : martensite ration in the beginning of loading cycle
        % zT0 : temprature induced martensite ration in the beginning of loading cycle
        % zS0 : Stress induced martensite ration in the beginning of loading cycle
        % Note that Ms < T
        % Note that this relation stands true only if: 
        % criticalStressStart + CM*(T-Ms) < stress < criticalStressFinisg + CM*(T-Ms)
        % Also note that phase transformation from austenite to detwinned martensite
        % happens only when specimen is in the "Loading" phase and stress is
        % increasing (0 < zDot)
    
        % Update: in this update we removed the criticalStressStart and
        % criticalStressFinish from the required arguments. the user has to
        % only provide S_M_start and S_M_finish (Critical stress that
        % conversion to austenite starts and ends respectivelt). the user has
        % to calculate these parameters and pass them to method. In some papers
        % there is no criticalStress provided (unlike the brinson's model. we use
        % this method to make the function more friendly for implementing other
        % papaers)
    
        % Update: The function will nolonger update z when we are forcing the
        % spring to act linear
        if (Ms < T) && (S_M_start < stress && stress < S_M_finish) && ~params.forceLinear
            zS = (1-zS0)/2*cos(pi/(S_M_start - S_M_finish) * (stress - S_M_finish)) + (1 + zS0)/2;
            zT = zT0 - zT0 / (1-zT0) * (zS - zS0);
            z = zS + zT;
        else
            zS = zS0;
            zT = zT0;
            z = zS + zT;
        end
    end

    function [z] = get_z(stress,stressDot,T,Ms,As,zS0,zT0,currZ, S_M_finish, S_M_start, S_A_finish, S_A_start, params)
        % This function returns total martensite ratio
        z = -99;
        if 0 < stressDot
            % Loading
            [z,~,~] = Func_Brinson_Conversion_To_Detwinned_Martensite(T,stress,Ms,zS0,zT0, S_M_finish, S_M_start, S_A_finish, S_A_start, params);
        elseif stressDot == 0
            z = -9999999;
        elseif stressDot < 0
            % Unloading
            [z,~,~] = Func_Brinson_Conversion_To_Austenite(T,stress,As,zS0,zT0,currZ, S_M_finish, S_M_start, S_A_finish, S_A_start, params);
        end
    end

    function stress = stressFromStrainArray(strain, straindDot, strain_A2M, strain_M2A, q_A2M, q_M2A)
        % This function gets thestress for each strain with respect to
        % direction of the strain. This function is used instead of the vpasolve
        % in the previouse version. This causes the speed of the application to
        % increas roughly 20 times.
    
        % Args: 
        % strain: double: -
        % strainDot: +1 or -1: direction of the loading, if +1 austenite is
        %   changing to martensinte and vice versa
        % strain_A2M, strain_M2A, q_A2M, q_M2A: matrix: the stress and strains
        %   that the transformation happens
    
        % Returns: double: The stress for the strain.
    
        % ** Should run when we need to find stress when we have a change from

        if 0 <= straindDot
            % Austenite to martensite
            upperBoundIDX = find(strain < strain_A2M(:),1,"last");
            lowerBoundIDX = find(strain_A2M(:) <= strain, 1,"first");

            if isempty(upperBoundIDX) || isempty(lowerBoundIDX)  
                stress = [];
            else
                stress = q_A2M(lowerBoundIDX) + (q_A2M(upperBoundIDX) - q_A2M(lowerBoundIDX)) / (strain_A2M(upperBoundIDX) - strain_A2M(lowerBoundIDX)) * (strain - strain_A2M(lowerBoundIDX));
            end
            
        else
            % Martensite to austenite
            upperBoundIDX = find(strain < strain_M2A(:),1,"last");
            lowerBoundIDX = find(strain_M2A(:) <= strain, 1,"first");
        
    
            if isempty(upperBoundIDX) || isempty(lowerBoundIDX)  
                stress = [];
            else
                stress = q_M2A(lowerBoundIDX) + (q_M2A(upperBoundIDX) - q_M2A(lowerBoundIDX)) / (strain_M2A(upperBoundIDX) - strain_M2A(lowerBoundIDX)) * (strain - strain_M2A(lowerBoundIDX));
            end
        
        end
    end
    
    function tbl =  strainTableGenerator(A2M_CORE, M2A_CORE, q_A2M, q_M2A, params)
        % This function calculates the respective strain for each particular
        % stress.
        % ** Should run every time z0 chanegs.
    
        % Args:
        % A2M_CORE, M2A_CORE: matrix: Matrices that contain the core of the
        %   calculations. cos(pi/(S_M_start - S_M_finish) * (stress - S_M_finish))
        %   for austenite to martensite transform and cos(pi/(S_M_start - S_M_finish) * (stress - S_M_finish))
        %   for martensite to austenite.
        % q_A2M, q_M2A: double: The stresses at which austenite to martensite phase transform happens
        % params: struct
    
        % Returns:
        % A table containig the respective strain for each stress in q_A2M or
        % q_M2A matrix.
        Z_A2M = (1-params.zS0)/2*A2M_CORE+(1+params.zS0)/2;
        strain_A2M = ((Z_A2M/params.Gm + (1-Z_A2M)/params.Ga)) .* q_A2M + params.eL * Z_A2M;
        
        Z_M2A = params.z0/2*M2A_CORE;
        strain_M2A = ((Z_M2A/params.Gm + (1-Z_M2A)/params.Ga)) .* q_M2A + params.eL * Z_M2A;
    
        tbl = table(strain_A2M, strain_M2A);
    end
    
    function [G] = Func_SMA_Shear_Modulus(z,Gm,Ga, method)
        %FUNC_SMA_SHEAR_MODULUS 
        % in this function we get martensite fraction, (Percentage of martensite in alloy) to
        % outPut SMA's Shear's modulus. 
        % option 1 => We use Reuss scheme (Auricchio and Sacco, 1999): G = (z/Gm + (1-z)/Ga)^-1
        % option 2 => We use Linear destribution: G = Ga+ z*(Gm - Ga)
        % Gm: martensit's Young's modulus
        % Ga: austenite's Young's modulus
        % ** Note that 0 < z < 1 **
    
    
        %**NOTE** using method 2 for calculation, has some convergance issues.
        %use method 1 in the code
         
        if (method == 1)
            G = (z/Gm+(1-z)/Ga)^(-1);
        elseif (method == 2)
            G = Ga+ z*(Gm - Ga);
        else
            disp("ERROR in calculating shear modulus, wrong option entered!")
            G = -1;
        end
        
    end
    
    
end