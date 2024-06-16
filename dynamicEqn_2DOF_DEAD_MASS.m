function [dydt, params] = dynamicEqn_2DOF_DEAD_MASS(t,Y,params)
    % Params has to be a struct, containing the names of the parameters and
    % their values. Also remember that because of the way the dynamic
    % equation is calculated, md has to be waight of the absorber, not the
    % mass. In calculating the lagrangian, we have used Weight of md
    % instead of its mass. This choice has been taken for performance
    % reasons.
    % Y = [th2, Dth2, th1, Dth1]
    params.md = params.md*params.g;
    m1 = params.m1;
    m2 = params.m2;
    g = params.g;
    md = params.md;
    params.md = md;
    c1 = params.c1;
    c2 = params.c2;
    kt1Mul = params.kt1Mul;
    kt2Mul = params.kt2Mul;
    l1 = params.l1;
    lc1 = params.lc1;
    lc2 = params.lc2;
    ld = params.ld;
    i1 = params.i1;
    i2 = params.i2;

    ampElbow = params.ampElbow;
    ampWrist = params.ampWrist ;
    freqElbow = params.freqElbow;
    freqWrist = params.freqWrist;


    [kt1,kt2,~] = jointStiffness(Y(1), 0, Y(4));
    kt1 = kt1 * kt1Mul;
    kt2 = kt2 * kt2Mul;
    
    excitationElbow = ampElbow * cos(t * freqElbow);
    excitationWrist = ampWrist * cos(t * freqWrist);
    [elbowTorque, wristTorque] = postureTorque(Y(3), Y(1), params); % The torques required to hold elbow and wrist still
    T1 = elbowTorque + excitationElbow;
    T2 = wristTorque + excitationWrist;

    dydt = [Y(2);-(T1*i2-T2*i1-T2*l1^2*m2-T2*lc1^2*m1+T1*lc2^2*m2-T2*ld^2*md+c2*i1*Y(2)+c2*i2*Y(2)-c1*i2*Y(4)+i1*kt2*Y(1)+i2*kt2*Y(1)-i2*kt1*Y(3)+c2*l1^2*m2*Y(2)+c2*lc1^2*m1*Y(2)+c2*lc2^2*m2*Y(2)-c1*lc2^2*m2*Y(4)+c2*ld^2*md*Y(2)+kt2*l1^2*m2*Y(1)+kt2*lc1^2*m1*Y(1)+kt2*lc2^2*m2*Y(1)-kt1*lc2^2*m2*Y(3)+kt2*ld^2*md*Y(1)+sin(Y(1))*l1*lc2^3*m2^2*Y(2)^2+sin(Y(1))*l1*lc2^3*m2^2*Y(4)^2+sin(Y(1))*l1^3*lc2*m2^2*Y(4)^2+cos(Y(1)+Y(3))*g*i1*lc2*m2-cos(Y(3))*g*l1*lc2^2*m2^2+cos(Y(1)+Y(3))*g*l1^2*lc2*m2^2+cos(Y(1))*T1*l1*lc2*m2-cos(Y(1))*T2*l1*lc2*m2-cos(Y(3))*g*i2*l1*m2-cos(Y(3))*g*i2*lc1*m1-cos(Y(3))*g*i2*ld*md+sin(Y(1))*i2*l1*lc2*m2*Y(2)^2+sin(Y(1))*i1*l1*lc2*m2*Y(4)^2+sin(Y(1))*i2*l1*lc2*m2*Y(4)^2-cos(Y(3))*g*lc1*lc2^2*m1*m2-cos(Y(3))*g*lc2^2*ld*m2*md+cos(Y(1)+Y(3))*g*lc1^2*lc2*m1*m2+cos(Y(1)+Y(3))*g*lc2*ld^2*m2*md+2*sin(Y(1))*l1*lc2^3*m2^2*Y(2)*Y(4)+cos(Y(1))*sin(Y(1))*l1^2*lc2^2*m2^2*Y(2)^2+2*cos(Y(1))*sin(Y(1))*l1^2*lc2^2*m2^2*Y(4)^2+2*cos(Y(1))*c2*l1*lc2*m2*Y(2)-cos(Y(1))*c1*l1*lc2*m2*Y(4)+2*cos(Y(1))*kt2*l1*lc2*m2*Y(1)-cos(Y(1))*kt1*l1*lc2*m2*Y(3)-cos(Y(1))*cos(Y(3))*g*l1^2*lc2*m2^2+cos(Y(1))*cos(Y(1)+Y(3))*g*l1*lc2^2*m2^2+sin(Y(1))*l1*lc1^2*lc2*m1*m2*Y(4)^2+sin(Y(1))*l1*lc2*ld^2*m2*md*Y(4)^2+2*sin(Y(1))*i2*l1*lc2*m2*Y(2)*Y(4)+2*cos(Y(1))*sin(Y(1))*l1^2*lc2^2*m2^2*Y(2)*Y(4)-cos(Y(1))*cos(Y(3))*g*l1*lc1*lc2*m1*m2-cos(Y(1))*cos(Y(3))*g*l1*lc2*ld*m2*md)/(i1*i2+l1^2*lc2^2*m2^2+i2*l1^2*m2+i2*lc1^2*m1+i1*lc2^2*m2+i2*ld^2*md-cos(Y(1))^2*l1^2*lc2^2*m2^2+lc1^2*lc2^2*m1*m2+lc2^2*ld^2*m2*md);Y(4);(T1*i2+T1*lc2^2*m2+c2*i2*Y(2)-c1*i2*Y(4)+i2*kt2*Y(1)-i2*kt1*Y(3)+c2*lc2^2*m2*Y(2)-c1*lc2^2*m2*Y(4)+kt2*lc2^2*m2*Y(1)-kt1*lc2^2*m2*Y(3)+sin(Y(1))*l1*lc2^3*m2^2*Y(2)^2+sin(Y(1))*l1*lc2^3*m2^2*Y(4)^2-cos(Y(3))*g*l1*lc2^2*m2^2-cos(Y(1))*T2*l1*lc2*m2-cos(Y(3))*g*i2*l1*m2-cos(Y(3))*g*i2*lc1*m1-cos(Y(3))*g*i2*ld*md+sin(Y(1))*i2*l1*lc2*m2*Y(2)^2+sin(Y(1))*i2*l1*lc2*m2*Y(4)^2-cos(Y(3))*g*lc1*lc2^2*m1*m2-cos(Y(3))*g*lc2^2*ld*m2*md+2*sin(Y(1))*l1*lc2^3*m2^2*Y(2)*Y(4)+cos(Y(1))*sin(Y(1))*l1^2*lc2^2*m2^2*Y(4)^2+cos(Y(1))*c2*l1*lc2*m2*Y(2)+cos(Y(1))*kt2*l1*lc2*m2*Y(1)+cos(Y(1))*cos(Y(1)+Y(3))*g*l1*lc2^2*m2^2+2*sin(Y(1))*i2*l1*lc2*m2*Y(2)*Y(4))/(i1*i2+l1^2*lc2^2*m2^2+i2*l1^2*m2+i2*lc1^2*m1+i1*lc2^2*m2+i2*ld^2*md-cos(Y(1))^2*l1^2*lc2^2*m2^2+lc1^2*lc2^2*m1*m2+lc2^2*ld^2*m2*md);];

    function [T_elbow, T_wrist] = postureTorque(th1, th2, Params)
        % For the wrist and hand to keep its initial posture, the body needs
        % to excert a certain amount of torque. This function calculates the
        % torques that will keep the wrist and elbow from falling.
    
        T_wrist = +Params.lc2*Params.g*Params.m2*cos(th2+th1); % same as T2
        T_elbow = Params.g*Params.m1*Params.lc1*cos(th1) + Params.m2*Params.g*(Params.l1*cos(th1)+Params.lc2*cos(th2+th1)) + Params.md*Params.g*(Params.ld*cos(th1)) - T_wrist; % same as T1
    end
    
    function [kt1,kt2,kt3] = jointStiffness(th2, th3, th1Dir)
        % Part1: This part gets the stifness of the elbow with respect to angle
        
        % Args:
        % th1_dir : float: The direction of elbows movemebt (Flexion or Extension)
        
        % Returns:
        % The torque antagonist torque generated in elbow in (N.m)
        
        th1Dir = sign(th1Dir);
    
        if 0 <= th1Dir
            kt1 = 0.0516*180/pi;
        else
            kt1 = 0.0516*180/pi;
        end
    
            % Part2: This part gets the stifness of the wrist with respect to two angles,
        % 1- FE (Flexion Extension)(th2) and 2- RUD (Radial Ulnar  Deviation)(th3). 
        % The reference for the Data in this function is : DOI 10.1152/jn.01014.2011
        % In this paper the stifness of the wrist for small angles (-15 to +15
        % degrees) has been obtained. with multiiplying the wrist stiffness to
        % the wrist angle, we can get the antagonist torque in the wrist.
    
        % We have 24 experimental datas, with an elliptic regression we have
        % acquired the formulation the ellipce that fits onto the experimental data.
        % The ellipse's formulation is as follows:
        % 0.4352x^2 + 0.2482xy + 0.2231y^2 + -0.1344x + -0.2594y = 1
        % To get the python code for fitting the ellipse, please see fitToEllipse.ipynb
    
        % To get the palm's we use the formula theta = atan(tan(alpha) / tan(beta)) 
        % which alpha is for FE and beta is for RUD.
    
        % Args:
        % FE: float: the flexion extension angle (alpha) in radians
        % RUD: float: the radial ulnar deviation angle (beta) in radians
    
        % Returns: 
        % k: float: The stiffness of the joint in the palms plain of motion (N.m/rad)
    
        alpha = th2; % FE
        beta = th3; % RUD
    
        theta = atan2(tan(alpha),tan(beta));
    
        % Here we have the angle theta, now we have to see where the line
        % passing from origin, and having the ngle theta with x axis,
        % intersects with regrtessed ellipse.
    
        R = tan(theta); % Hence x = R * y (refere to the chart)
    
        % The equation is as follows. we have to solve it: 
        % eqn = 0.4352*(R*y)^2 + 0.2482*(R*y)*y + 0.2231*y^2 + -0.1344*R*y + -0.2594*y == 1;
        coeff1 = 0.4352*R^2 +  0.2482*R + 0.2231; % Coefficient for y^2
        coeff2 = -0.1344*R -0.2594; % Coefficient for y^1
        coeff3 = -1; % Coefficient for y^0
    
        sol = roots([coeff1, coeff2, coeff3]);
    
        % We always have two solutions. if beta is positive, then we take the
        % positive answer as solution. if not, the negative solution. Also note
        % that always, one solution is positive and another negative.
        if 0 <= beta
            y = max(real(sol));
        else
            y = min(real(sol));
        end
    
        if alpha == 0  && beta == 0
            k2 = 0;
            %fprintf("k=%f\n", k);
        else
            x = R * y;
            k2 = sqrt(x^2+y^2);
            %fprintf("alpha =%f \nbeta=%f \ntheta=%f   \nx=%f   \ny=%f  \nk=%f\n------------\n",alpha,beta, theta*180/pi, x, y, k);
        end
        
        phi = atan(sqrt(tan(alpha)^2+tan(beta)^2));
        %fprintf("alpha =%f \nbeta=%f \ntheta=%f   \nphi=%f   \nk=%f\n------------\n",alpha,beta, theta*180/pi,phi*180/pi, k);
        Torque = k2 * phi;
        Torque_FE  = Torque * sin(theta);
        Torque_RUD = Torque * cos(theta);
        if th2 == 0; kt2 = 0; else ;kt2 = Torque_FE / th2; end
        if th3 == 0; kt3 = 0; else ;kt3 = Torque_RUD / th3; end

    end
    params.md = params.md/params.g; % Because we pass params back to ode15s solver, we have to change md to kg so we avoid accumulation of gravity multipliers in absorber mass
end
