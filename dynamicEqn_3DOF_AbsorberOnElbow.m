function [dydt, params] = dynamicEqn_3DOF_AbsorberOnElbow(t,Y,params)
    % Params has to be a struct, containing the names of the parameters and
    % their values. Also remember that because of the way the dynamic
    % equation is calculated, md has to be waight of the absorber, not the
    % mass. In calculating the lagrangian, we have used Weight of md
    % instead of its mass. This choice has been taken for performance
    % reasons.
    % Y = [th2, Dth2, th1, Dth1, xd, Dxd]
    params.md = params.md*params.g;
%     m1 = params.m1;
%     m2 = params.m2;
%     c1 = params.c1;
%     c2 = params.c2;
%     g = params.g;
    md = params.md;
    kt1Mul = params.kt1Mul;
    kt2Mul = params.kt2Mul;
    cd = params.cd;
    kd = params.kd;
    ld = params.ld;
%     l1 = params.l1;
%     lc1 = params.lc1;
%     lc2 = params.lc2;
%     i1 = params.i1;
%     i2 = params.i2;

    ampElbow = params.ampElbow;
    ampWrist = params.ampWrist;
    freqElbow = params.freqElbow;
    freqWrist= params.freqWrist;
    
    [kt1,kt2,~] = jointStiffness(Y(1), 0, Y(4));
    kt1 = kt1 * kt1Mul;
    kt2 = kt2 * kt2Mul;
    
    excitationElbow = ampElbow * cos(t * freqElbow);
    excitationWrist = ampWrist * cos(t * freqWrist);

    [elbowTorque, wristTorque] = postureTorque(Y(3), Y(1), Y(5), params); % The torques required to hold elbow and wrist still
    params.md = params.md/params.g; % Because we pass params back to ode15s solver, we have to change md to kg so we avoid accumulation of gravity multipliers in absorber mass
    T1 = elbowTorque + excitationElbow;
    T2 = wristTorque + excitationWrist;

    F_SMA = 0;
    if params.spring1_Active || params.spring2_Active
        [params.spring1,~,spring1_force] = SMA_Spring_Excite(params.spring1_temp, (params.spring1_temp.preStretchDisplacement + Y(5)) * 1000);
        [params.spring2,~,spring2_force] = SMA_Spring_Excite(params.spring2_temp, (params.spring2_temp.preStretchDisplacement - Y(5)) * 1000); 
        spring2_force = -1 * spring2_force;

        totlForce  = (spring1_force + spring2_force);
        % I dont know why but when i use F_SMA, the amplitude when the
        % frequency of forces are zero differs from the system of simple
        % mass and damper (Even when there is no transformation and the 
        % sma spring is essentially a linear spring made of asutenite).
        % Even when i replace the F_SMA with kd*x and assign 0 to kd, that 
        % should result in a system with linear spring, i get different 
        % results than a system with a linear spring and no F_SMA. I think
        % this is due to SMA's force being conservative so i calculate the 
        % approximate kd every time i want to use F_SMA(but SMA's force is
        % path dependant so it shouldnt be conservative)
        if Y(5) == 0
            % When we have no displacement, the spring is full austenite.
            % we do this to avoid division by zero
            kd = params.kd;% + spring1_temp.Ka;
        else 
            kd = params.kd + abs(totlForce/Y(5));
        end
%         F_SMA = totlForce;

    end
    
    % The main equation is below but for performance reasons, we have done
    % the redundant multiplications; The second equation matrix is only a
    % function of T1 T2 F_SMA k1 k2 and other design variables that change
    % during the optimization.
%     dydt = [Y(2);-(T1*i2-T2*i1-T2*md*Y(5)^2-T2*l1^2*m2-T2*lc1^2*m1+T1*lc2^2*m2+c2*i1*Y(2)+c2*i2*Y(2)-c1*i2*Y(4)+i1*kt2*Y(1)+i2*kt2*Y(1)-i2*kt1*Y(3)+c2*l1^2*m2*Y(2)+c2*lc1^2*m1*Y(2)+c2*lc2^2*m2*Y(2)-c1*lc2^2*m2*Y(4)+kt2*l1^2*m2*Y(1)+kt2*lc1^2*m1*Y(1)+kt2*lc2^2*m2*Y(1)-kt1*lc2^2*m2*Y(3)+cd*i2*ld*Y(6)+i2*kd*ld*Y(5)+c2*md*Y(2)*Y(5)^2+kt2*md*Y(1)*Y(5)^2+sin(Y(1))*l1*lc2^3*m2^2*Y(2)^2+sin(Y(1))*l1*lc2^3*m2^2*Y(4)^2+sin(Y(1))*l1^3*lc2*m2^2*Y(4)^2+cos(Y(1)+Y(3))*g*i1*lc2*m2-2*i2*md*Y(4)*Y(5)*Y(6)-cos(Y(3))*g*l1*lc2^2*m2^2+cos(Y(1)+Y(3))*g*l1^2*lc2*m2^2-i2*ld*md*Y(4)^2*Y(5)+sin(Y(3))*g*i2*md*Y(5)+cos(Y(1))*T1*l1*lc2*m2-cos(Y(1))*T2*l1*lc2*m2+cd*lc2^2*ld*m2*Y(6)-cos(Y(3))*g*i2*l1*m2-cos(Y(3))*g*i2*lc1*m1+kd*lc2^2*ld*m2*Y(5)-lc2^2*ld*m2*md*Y(4)^2*Y(5)+sin(Y(3))*g*lc2^2*m2*md*Y(5)+sin(Y(1))*i2*l1*lc2*m2*Y(2)^2+sin(Y(1))*i1*l1*lc2*m2*Y(4)^2+sin(Y(1))*i2*l1*lc2*m2*Y(4)^2+cos(Y(1)+Y(3))*g*lc2*m2*md*Y(5)^2-cos(Y(3))*g*lc1*lc2^2*m1*m2+cos(Y(1)+Y(3))*g*lc1^2*lc2*m1*m2+2*sin(Y(1))*l1*lc2^3*m2^2*Y(2)*Y(4)+cos(Y(1))*sin(Y(1))*l1^2*lc2^2*m2^2*Y(2)^2+2*cos(Y(1))*sin(Y(1))*l1^2*lc2^2*m2^2*Y(4)^2-2*lc2^2*m2*md*Y(4)*Y(5)*Y(6)+2*cos(Y(1))*c2*l1*lc2*m2*Y(2)-cos(Y(1))*c1*l1*lc2*m2*Y(4)+2*cos(Y(1))*kt2*l1*lc2*m2*Y(1)-cos(Y(1))*kt1*l1*lc2*m2*Y(3)-cos(Y(1))*cos(Y(3))*g*l1^2*lc2*m2^2+cos(Y(1))*cos(Y(1)+Y(3))*g*l1*lc2^2*m2^2+sin(Y(1))*l1*lc2*m2*md*Y(4)^2*Y(5)^2+sin(Y(1))*l1*lc1^2*lc2*m1*m2*Y(4)^2+2*sin(Y(1))*i2*l1*lc2*m2*Y(2)*Y(4)+cos(Y(1))*cd*l1*lc2*ld*m2*Y(6)+cos(Y(1))*kd*l1*lc2*ld*m2*Y(5)+2*cos(Y(1))*sin(Y(1))*l1^2*lc2^2*m2^2*Y(2)*Y(4)-2*cos(Y(1))*l1*lc2*m2*md*Y(4)*Y(5)*Y(6)-cos(Y(1))*l1*lc2*ld*m2*md*Y(4)^2*Y(5)+cos(Y(1))*sin(Y(3))*g*l1*lc2*m2*md*Y(5)-cos(Y(1))*cos(Y(3))*g*l1*lc1*lc2*m1*m2)/(i1*i2+l1^2*lc2^2*m2^2+i2*md*Y(5)^2+i2*l1^2*m2+i2*lc1^2*m1+i1*lc2^2*m2+lc2^2*m2*md*Y(5)^2-cos(Y(1))^2*l1^2*lc2^2*m2^2+lc1^2*lc2^2*m1*m2);Y(4);(T1*i2+T1*lc2^2*m2+c2*i2*Y(2)-c1*i2*Y(4)+i2*kt2*Y(1)-i2*kt1*Y(3)+c2*lc2^2*m2*Y(2)-c1*lc2^2*m2*Y(4)+kt2*lc2^2*m2*Y(1)-kt1*lc2^2*m2*Y(3)+cd*i2*ld*Y(6)+i2*kd*ld*Y(5)+sin(Y(1))*l1*lc2^3*m2^2*Y(2)^2+sin(Y(1))*l1*lc2^3*m2^2*Y(4)^2-2*i2*md*Y(4)*Y(5)*Y(6)-cos(Y(3))*g*l1*lc2^2*m2^2-i2*ld*md*Y(4)^2*Y(5)+sin(Y(3))*g*i2*md*Y(5)-cos(Y(1))*T2*l1*lc2*m2+cd*lc2^2*ld*m2*Y(6)-cos(Y(3))*g*i2*l1*m2-cos(Y(3))*g*i2*lc1*m1+kd*lc2^2*ld*m2*Y(5)-lc2^2*ld*m2*md*Y(4)^2*Y(5)+sin(Y(3))*g*lc2^2*m2*md*Y(5)+sin(Y(1))*i2*l1*lc2*m2*Y(2)^2+sin(Y(1))*i2*l1*lc2*m2*Y(4)^2-cos(Y(3))*g*lc1*lc2^2*m1*m2+2*sin(Y(1))*l1*lc2^3*m2^2*Y(2)*Y(4)+cos(Y(1))*sin(Y(1))*l1^2*lc2^2*m2^2*Y(4)^2-2*lc2^2*m2*md*Y(4)*Y(5)*Y(6)+cos(Y(1))*c2*l1*lc2*m2*Y(2)+cos(Y(1))*kt2*l1*lc2*m2*Y(1)+cos(Y(1))*cos(Y(1)+Y(3))*g*l1*lc2^2*m2^2+2*sin(Y(1))*i2*l1*lc2*m2*Y(2)*Y(4))/(i1*i2+l1^2*lc2^2*m2^2+i2*md*Y(5)^2+i2*l1^2*m2+i2*lc1^2*m1+i1*lc2^2*m2+lc2^2*m2*md*Y(5)^2-cos(Y(1))^2*l1^2*lc2^2*m2^2+lc1^2*lc2^2*m1*m2);Y(6);-(i2*kd*md*Y(5)^3-F_SMA*i1*i2-F_SMA*i2*l1^2*m2-F_SMA*i2*lc1^2*m1-F_SMA*i1*lc2^2*m2+cd*i1*i2*Y(6)-i2*md^2*Y(4)^2*Y(5)^3+i1*i2*kd*Y(5)+T1*i2*ld*md-F_SMA*l1^2*lc2^2*m2^2-F_SMA*i2*md*Y(5)^2+T1*lc2^2*ld*m2*md-F_SMA*lc2^2*m2*md*Y(5)^2+cos(Y(1))^2*F_SMA*l1^2*lc2^2*m2^2+kd*lc2^2*m2*md*Y(5)^3-F_SMA*lc1^2*lc2^2*m1*m2+cos(Y(3))*g*i2*md^2*Y(5)^2-i2*ld^2*md^2*Y(4)^2*Y(5)+c2*i2*ld*md*Y(2)-c1*i2*ld*md*Y(4)+i2*kt2*ld*md*Y(1)-i2*kt1*ld*md*Y(3)+cd*l1^2*lc2^2*m2^2*Y(6)+kd*l1^2*lc2^2*m2^2*Y(5)+cd*i2*md*Y(5)^2*Y(6)-i1*i2*md*Y(4)^2*Y(5)+cd*i2*l1^2*m2*Y(6)+cd*i2*lc1^2*m1*Y(6)+cd*i1*lc2^2*m2*Y(6)+cd*i2*ld^2*md*Y(6)-lc2^2*m2*md^2*Y(4)^2*Y(5)^3+cos(Y(3))*g*i1*i2*md+i2*kd*l1^2*m2*Y(5)+i2*kd*lc1^2*m1*Y(5)+i1*kd*lc2^2*m2*Y(5)+i2*kd*ld^2*md*Y(5)+cd*lc2^2*m2*md*Y(5)^2*Y(6)-i2*l1^2*m2*md*Y(4)^2*Y(5)-i2*lc1^2*m1*md*Y(4)^2*Y(5)-i1*lc2^2*m2*md*Y(4)^2*Y(5)+sin(Y(3))*g*i2*ld*md^2*Y(5)-cos(Y(1))^2*cd*l1^2*lc2^2*m2^2*Y(6)+cd*lc1^2*lc2^2*m1*m2*Y(6)+cd*lc2^2*ld^2*m2*md*Y(6)-cos(Y(1))^2*kd*l1^2*lc2^2*m2^2*Y(5)+cos(Y(3))*g*i2*l1^2*m2*md+cos(Y(3))*g*i2*lc1^2*m1*md+cos(Y(3))*g*i1*lc2^2*m2*md+kd*lc1^2*lc2^2*m1*m2*Y(5)+kd*lc2^2*ld^2*m2*md*Y(5)-2*i2*ld*md^2*Y(4)*Y(5)*Y(6)+cos(Y(3))*g*lc2^2*m2*md^2*Y(5)^2-l1^2*lc2^2*m2^2*md*Y(4)^2*Y(5)-lc2^2*ld^2*m2*md^2*Y(4)^2*Y(5)+c2*lc2^2*ld*m2*md*Y(2)-c1*lc2^2*ld*m2*md*Y(4)+kt2*lc2^2*ld*m2*md*Y(1)-kt1*lc2^2*ld*m2*md*Y(3)+cos(Y(3))*g*l1^2*lc2^2*m2^2*md+cos(Y(1))^2*l1^2*lc2^2*m2^2*md*Y(4)^2*Y(5)-lc1^2*lc2^2*m1*m2*md*Y(4)^2*Y(5)+sin(Y(3))*g*lc2^2*ld*m2*md^2*Y(5)-cos(Y(1))^2*cos(Y(3))*g*l1^2*lc2^2*m2^2*md-cos(Y(3))*g*l1*lc2^2*ld*m2^2*md+cos(Y(3))*g*lc1^2*lc2^2*m1*m2*md-cos(Y(1))*T2*l1*lc2*ld*m2*md-cos(Y(3))*g*i2*l1*ld*m2*md-cos(Y(3))*g*i2*lc1*ld*m1*md+sin(Y(1))*l1*lc2^3*ld*m2^2*md*Y(2)^2+sin(Y(1))*l1*lc2^3*ld*m2^2*md*Y(4)^2-2*lc2^2*ld*m2*md^2*Y(4)*Y(5)*Y(6)+2*sin(Y(1))*l1*lc2^3*ld*m2^2*md*Y(2)*Y(4)+cos(Y(1))*sin(Y(1))*l1^2*lc2^2*ld*m2^2*md*Y(4)^2+cos(Y(1))*c2*l1*lc2*ld*m2*md*Y(2)+cos(Y(1))*kt2*l1*lc2*ld*m2*md*Y(1)+cos(Y(1))*cos(Y(1)+Y(3))*g*l1*lc2^2*ld*m2^2*md+sin(Y(1))*i2*l1*lc2*ld*m2*md*Y(2)^2+sin(Y(1))*i2*l1*lc2*ld*m2*md*Y(4)^2-cos(Y(3))*g*lc1*lc2^2*ld*m1*m2*md+2*sin(Y(1))*i2*l1*lc2*ld*m2*md*Y(2)*Y(4))/(md*(i1*i2+l1^2*lc2^2*m2^2+i2*md*Y(5)^2+i2*l1^2*m2+i2*lc1^2*m1+i1*lc2^2*m2+lc2^2*m2*md*Y(5)^2-cos(Y(1))^2*l1^2*lc2^2*m2^2+lc1^2*lc2^2*m1*m2));];    
    dydt = [Y(2);-(1.0*(0.0010713042*Y(2) - 0.00034965485*Y(4) + 0.0034965485*T1 - 0.050068663*T2 + 0.015351102*cos(Y(3) + Y(1)) - 0.008571361*cos(Y(3)) - 0.01925997*cos(Y(3))*cos(Y(1)) - 0.0034965485*kt1*Y(3) + 0.053565211*kt2*Y(1) - 0.0007856794*Y(4)*cos(Y(1)) + 0.00031427176*Y(2)*cos(Y(1)) + 0.007856794*T1*cos(Y(1)) - 0.007856794*T2*cos(Y(1)) + 0.0024089008*cos(Y(3) + Y(1))*cos(Y(1)) + 0.00042085083*Y(4)^2*sin(Y(1)) + 0.000027471661*Y(2)^2*sin(Y(1)) - 0.007856794*kt1*Y(3)*cos(Y(1)) + 0.015713588*kt2*Y(1)*cos(Y(1)) + 0.034289427*md*Y(5)*sin(Y(3)) + 0.30660099*md*Y(5)^2*cos(Y(3) + Y(1)) + 0.0034965485*Y(6)*cd*ld + 0.0034965485*kd*ld*Y(5) + 0.00012345842*Y(4)^2*cos(Y(1))*sin(Y(1)) + 0.000061729212*Y(2)^2*cos(Y(1))*sin(Y(1)) + 0.000054943322*Y(4)*Y(2)*sin(Y(1)) + 0.02*Y(2)*md*Y(5)^2 - 1.0*T2*md*Y(5)^2 + 0.077048829*md*Y(5)*cos(Y(1))*sin(Y(3)) + 0.007856794*Y(4)^2*md*Y(5)^2*sin(Y(1)) + 0.007856794*Y(6)*cd*ld*cos(Y(1)) - 0.0034965485*Y(4)^2*ld*md*Y(5) + kt2*md*Y(1)*Y(5)^2 + 0.007856794*kd*ld*Y(5)*cos(Y(1)) + 0.00012345842*Y(4)*Y(2)*cos(Y(1))*sin(Y(1)) - 0.0069930969*Y(4)*Y(6)*md*Y(5) - 0.007856794*Y(4)^2*ld*md*Y(5)*cos(Y(1)) - 0.015713588*Y(4)*Y(6)*md*Y(5)*cos(Y(1))))/(0.0034965485*md*Y(5)^2 - 0.000061729212*cos(Y(1))^2 + 0.00017506751);Y(4);(0.000069930969*Y(2) - 0.00034965485*Y(4) + 0.0034965485*T1 - 0.008571361*cos(Y(3)) - 0.0034965485*kt1*Y(3) + 0.0034965485*kt2*Y(1) + 0.00015713588*Y(2)*cos(Y(1)) - 0.007856794*T2*cos(Y(1)) + 0.0024089008*cos(Y(3) + Y(1))*cos(Y(1)) + 0.000027471661*Y(4)^2*sin(Y(1)) + 0.000027471661*Y(2)^2*sin(Y(1)) + 0.007856794*kt2*Y(1)*cos(Y(1)) + 0.034289427*md*Y(5)*sin(Y(3)) + 0.0034965485*Y(6)*cd*ld + 0.0034965485*kd*ld*Y(5) + 0.000061729212*Y(4)^2*cos(Y(1))*sin(Y(1)) + 0.000054943322*Y(4)*Y(2)*sin(Y(1)) - 0.0034965485*Y(4)^2*ld*md*Y(5) - 0.0069930969*Y(4)*Y(6)*md*Y(5))/(0.0034965485*md*Y(5)^2 - 0.000061729212*cos(Y(1))^2 + 0.00017506751);Y(6);-(1.0*(0.00017506751*Y(6)*cd - 0.00017506751*F_SMA + 0.00017506751*kd*Y(5) + 0.0017168258*md*cos(Y(3)) + 0.000061729212*F_SMA*cos(Y(1))^2 + 0.0034965485*kd*md*Y(5)^3 - 0.008571361*ld*md*cos(Y(3)) - 0.000061729212*Y(6)*cd*cos(Y(1))^2 - 0.000061729212*kd*Y(5)*cos(Y(1))^2 - 0.00034965485*Y(4)*ld*md + 0.000069930969*Y(2)*ld*md + 0.0034965485*T1*ld*md - 0.0034965485*Y(4)^2*md^2*Y(5)^3 + 0.034289427*md^2*Y(5)^2*cos(Y(3)) - 0.00060535677*md*cos(Y(3))*cos(Y(1))^2 - 0.00017506751*Y(4)^2*md*Y(5) - 0.0034965485*F_SMA*md*Y(5)^2 - 0.0034965485*Y(4)^2*ld^2*md^2*Y(5) - 0.0034965485*kt1*ld*md*Y(3) + 0.0034965485*kt2*ld*md*Y(1) + 0.000061729212*Y(4)^2*md*Y(5)*cos(Y(1))^2 + 0.0034965485*Y(6)*cd*ld^2*md + 0.0034965485*Y(6)*cd*md*Y(5)^2 + 0.00015713588*Y(2)*ld*md*cos(Y(1)) - 0.007856794*T2*ld*md*cos(Y(1)) + 0.0034965485*kd*ld^2*md*Y(5) + 0.0024089008*ld*md*cos(Y(3) + Y(1))*cos(Y(1)) + 0.000027471661*Y(4)^2*ld*md*sin(Y(1)) + 0.000027471661*Y(2)^2*ld*md*sin(Y(1)) + 0.034289427*ld*md^2*Y(5)*sin(Y(3)) - 0.0069930969*Y(4)*Y(6)*ld*md^2*Y(5) + 0.000061729212*Y(4)^2*ld*md*cos(Y(1))*sin(Y(1)) + 0.000054943322*Y(4)*Y(2)*ld*md*sin(Y(1)) + 0.007856794*kt2*ld*md*Y(1)*cos(Y(1))))/(md*(0.0034965485*md*Y(5)^2 - 0.000061729212*cos(Y(1))^2 + 0.00017506751))];
    
    function [T_elbow, T_wrist] = postureTorque(th1, th2, xd, Params)
        % For the wrist and hand to keep its initial posture, the body needs
        % to excert a certain amount of torque. This function calculates the
        % torques that will keep the wrist and elbow from falling.

        T_wrist = Params.lc2*Params.g*Params.m2*cos(th2+th1); % same as T2
        T_elbow = Params.g*Params.m1*Params.lc1*cos(th1) + Params.m2*Params.g*(Params.l1*cos(th1)+Params.lc2*cos(th2+th1)) + Params.md*Params.g*(Params.ld*cos(th1)+(xd+Params.xdInit)*cos(th1+pi/2)) - T_wrist; % same as T1 
    end
    
    function [kt1,kt2,kt3] = jointStiffness(th2, th3, th1Dir)
    %     % Part1: This part gets the stifness of the elbow with respect to angle
    %     
    %     % Args:
    %     % th1_dir : float: The direction of elbows movemebt (Flexion or Extension)
    %     
    %     % Returns:
    %     % The torque antagonist torque generated in elbow in (N.m)
    %     
    %     th1Dir = sign(th1Dir);
    % 
    %     if 0 <= th1Dir
    %         kt1 = 0.0516*180/pi;
    %     else
    %         kt1 = 0.0516*180/pi;
    %     end
    % 
    %         % Part2: This part gets the stifness of the wrist with respect to two angles,
    %     % 1- FE (Flexion Extension)(th2) and 2- RUD (Radial Ulnar  Deviation)(th3). 
    %     % The reference for the Data in this function is : DOI 10.1152/jn.01014.2011
    %     % In this paper the stifness of the wrist for small angles (-15 to +15
    %     % degrees) has been obtained. with multiiplying the wrist stiffness to
    %     % the wrist angle, we can get the antagonist torque in the wrist.
    % 
    %     % We have 24 experimental datas, with an elliptic regression we have
    %     % acquired the formulation the ellipce that fits onto the experimental data.
    %     % The ellipse's formulation is as follows:
    %     % 0.4352x^2 + 0.2482xy + 0.2231y^2 + -0.1344x + -0.2594y = 1
    %     % To get the python code for fitting the ellipse, please see fitToEllipse.ipynb
    % 
    %     % To get the palm's we use the formula theta = atan(tan(alpha) / tan(beta)) 
    %     % which alpha is for FE and beta is for RUD.
    % 
    %     % Args:
    %     % FE: float: the flexion extension angle (alpha) in radians
    %     % RUD: float: the radial ulnar deviation angle (beta) in radians
    % 
    %     % Returns: 
    %     % k: float: The stiffness of the joint in the palms plain of motion (N.m/rad)
    % 
    %     alpha = th2; % FE
    %     beta = th3; % RUD
    % 
    %     theta = atan2(tan(alpha),tan(beta));
    % 
    %     % Here we have the angle theta, now we have to see where the line
    %     % passing from origin, and having the ngle theta with x axis,
    %     % intersects with regrtessed ellipse.
    % 
    %     R = tan(theta); % Hence x = R * y (refere to the chart)
    % 
    %     % The equation is as follows. we have to solve it: 
    %     % eqn = 0.4352*(R*y)^2 + 0.2482*(R*y)*y + 0.2231*y^2 + -0.1344*R*y + -0.2594*y == 1;
    %     coeff1 = 0.4352*R^2 +  0.2482*R + 0.2231; % Coefficient for y^2
    %     coeff2 = -0.1344*R -0.2594; % Coefficient for y^1
    %     coeff3 = -1; % Coefficient for y^0
    % 
    %     sol = roots([coeff1, coeff2, coeff3]);
    % 
    %     % We always have two solutions. if beta is positive, then we take the
    %     % positive answer as solution. if not, the negative solution. Also note
    %     % that always, one solution is positive and another negative.
    %     if 0 <= beta
    %         y = max(real(sol));
    %     else
    %         y = min(real(sol));
    %     end
    % 
    %     if alpha == 0  && beta == 0
    %         k2 = 0;
    %         %fprintf("k=%f\n", k);
    %     else
    %         x = R * y;
    %         k2 = sqrt(x^2+y^2);
    %         %fprintf("alpha =%f \nbeta=%f \ntheta=%f   \nx=%f   \ny=%f  \nk=%f\n------------\n",alpha,beta, theta*180/pi, x, y, k);
    %     end
    %     
    %     phi = atan(sqrt(tan(alpha)^2+tan(beta)^2));
    %     %fprintf("alpha =%f \nbeta=%f \ntheta=%f   \nphi=%f   \nk=%f\n------------\n",alpha,beta, theta*180/pi,phi*180/pi, k);
    %     Torque = k2 * phi;
    %     Torque_FE  = Torque * sin(theta);
    %     Torque_RUD = Torque * cos(theta);
    %     if th2 == 0; kt2 = 0; else ;kt2 = Torque_FE / th2; end
    %     if th3 == 0; kt3 = 0; else ;kt3 = Torque_RUD / th3; end
        
        kt1 = 2.9564;
        kt2 = 1.6781;
        kt3 = 0;
    end
    

end