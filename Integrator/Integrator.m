classdef Integrator
% for integration of trajectories around small bodies. 
    
    properties
        SRPModel                % solar radiation pressure model
        thirdBodyGravityModel   % modelling background fields
        gravityModel            % bff gravity model of primary
        omega                   % rotation rate (rad/sec)
        
        integrator              % function handle for integrator
        mode                    % 0-inertial, frame 1-BFF, 2-variational Eqns
        odeOptions              % options from odeset
    end
    
    methods
        function obj = Integrator(gravityModel,integratorHandle)
            
            obj.mode = 0;
            if nargin == 1
                obj.integrator = @ode45;
                obj.gravityModel = gravityModel;
            elseif nargin == 2
                obj.gravityModel = gravityModel;
                obj.integrator = integratorHandle;
            elseif nargin~= 0 
                error('Incorrect number of inputs')
            end
            
        end
        function obj = setFrame(obj,frame)
            if  strcmp(frame,'BFF') || strcmp(frame,'bff')...
                    || strcmp(frame,'body-fixed-frame')
                obj.mode = 1;
            elseif strcmp(frame,'inertial')
                obj.mode = 0;
            else
                error("frame is 'inertial' or 'BFF'")
            end
        end
        function obj = setOdeOptions(obj,options)
            obj.odeOptions = options;
        end
                
        function [tout,xout] = integrate(obj,T,Xo)
            switch obj.mode
                case 0  
                    [tout,xout]=obj.integrator(@obj.GoverningEquations,T,Xo,obj.odeOptions);
                case 1  
                    [tout,xout]=obj.integrator(@obj.GoverningEquationsBFF,T,Xo,obj.odeOptions);
                case 2
                    [tout,xout]=obj.integrator(@obj.GoverningEquationsVariational,T,Xo,obj.odeOptions);
            end
        
        end
        function dx = GoverningEquations(obj, t, x )
        % governing equation for the inertial frame (default). 
        %------------------------------------------------------------------
        % Rotates pos to BFF calcs acceleration and rotates back
        %------------------------------------------------------------------
        % Inputs:
        %   t -- time
        %   x -- state vector
        %------------------------------------------------------------------
        % Outputs:
        %   dx - state vector derives
        %------------------------------------------------------------------
            
            %obj.ThirdBodyGravityModels([x(1),x(3),x(5)])

            Rot = [cos(obj.omega*t), -sin(obj.omega*t), 0;
                   sin(obj.omega*t),  cos(obj.omega*t), 0;
                   0,                 0,                1];
            
            P = (Rot'*[x(1),x(3),x(5)]')'; 
            
            [ a ] = obj.gravityModel.acceleration(P);
            
            a = (Rot*a')'; 
            
            dx = zeros(6,1);
            dx(1) = x(2);
            dx(2) = a(1);
            dx(3) = x(4);
            dx(4) = a(2);
            dx(5) = x(6);
            dx(6) = a(3);
            
        end
        function dx = GoverningEquationsBFF( obj, t, x )
        % governing equation for the body-fixed-frame. 
        %------------------------------------------------------------------
        % This will typically be slower than the default governing
        % equations because the orbital period is long and BFF results in
        % high curvature orbits compared to inertial frame
        %------------------------------------------------------------------
        % Inputs:
        %   t -- time
        %   x -- state vector
        %------------------------------------------------------------------
        % Outputs:
        %   dx - state vector derives
        %------------------------------------------------------------------
        
            [ a ] = obj.gravityModel.acceleration([x(1),x(3),x(5)]);
            
            dx = zeros(6,1);
            dx(1) = x(2);
            dx(2) = a(1)+obj.omega^2*x(1)+2*obj.omega*x(4);
            dx(3) = x(4);
            dx(4) = a(2)+obj.omega^2*x(3)-2*obj.omega*x(2);
            dx(5) = x(6);
            dx(6) = a(3);
            

        end
        function dx = GoverningEquationsVariational( obj, t, x)
            
            Rot = [cos(obj.omega*t), -sin(obj.omega*t), 0;
                   sin(obj.omega*t),  cos(obj.omega*t), 0;
                   0,             0,            1];
            
            P = (Rot'*[x(1),x(3),x(5)]')'; 
            
            a = obj.gravityModel.acceleration(P);
            gradA = obj.gravityModel.gravityGradient(P);
            
            gradA = (Rot*gradA*Rot');
            a = (Rot*a')';
            
            dx = zeros(12,1);
            dx(1) = x(2);
            dx(2) = a(1);
            dx(3) = x(4);
            dx(4) = a(2);
            dx(5) = x(6);
            dx(6) = a(3);
            dx(7) = x(8);
            dx(8) = gradA(1:3)*[x(7),x(9),x(11)]';
            dx(9) = x(10);
            dx(10) = gradA([2,4,5])*[x(7),x(9),x(11)]';
            dx(11) = x(12);
            dx(12) = gradA([3,5,6])*[x(7),x(9),x(11)]';
        end
        
    end
end

