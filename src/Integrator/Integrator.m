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
        % constructor
        %------------------------------------------------------------------
        % Inputs:
        %   gravityModel ------ gamma gravityModel object
        %   integratorhandle -- handle for matlab integrator func (default
        %                       ode45)
        %------------------------------------------------------------------   
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
        % method to swtich frame (body-fixed vs inrtial)
        %------------------------------------------------------------------
        % Inputs:
        %   frame -- string 'bff' or 'inertial'
        %------------------------------------------------------------------
            frame = lower(frame);
            if strcmp(frame,'bff')...
            || strcmp(frame,'body-fixed frame')
                obj.mode = 1;
            elseif strcmp(frame,'inertial')
                obj.mode = 0;
            else
                error("frame is 'inertial' or 'BFF'")
            end
        end
        function obj = setOdeOptions(obj,options)
        % pass-through for matlab standard integrator options
        %------------------------------------------------------------------
        % Inputs:
        %   options -- integrator options in matlab's default form
        %------------------------------------------------------------------
            obj.odeOptions = options;
        end
                
        function [tout,xout] = integrate(obj,T,Xo)
        % wrapper method 
        %------------------------------------------------------------------
        % Inputs:
        %   T -- [t0,t1] start and end time for integration
        %   Xo - state vector [x1,vx1,y1,vy1,z1,vz1,x2,vx2,...]
        %------------------------------------------------------------------
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

            % Rotate Position to body-fixed frame
            Rot = [cos(obj.omega*t), -sin(obj.omega*t), 0;
                   sin(obj.omega*t),  cos(obj.omega*t), 0;
                   0,                 0,                1];
            
            P = [x(1:6:end),x(3:6:end),x(5:6:end)]*Rot; 
            
            [ a ] = obj.gravityModel.acceleration(P);
            
            % rotate acceleration to inertial frame
            a = a*Rot'; 
            
            dx = zeros(size(x,1),1);
            dx(1:6:end) = x(2:6:end);
            dx(2:6:end) = a(:,1);
            dx(3:6:end) = x(4:6:end);
            dx(4:6:end) = a(:,2);
            dx(5:6:end) = x(6:6:end);
            dx(6:6:end) = a(:,3);
            
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
        
            [ a ] = obj.gravityModel.acceleration([x(1:6:end),x(3:6:end),x(5:6:end)]);
            
            rotTermx = obj.omega^2*x(1)+2*obj.omega*x(4);
            rotTermy = obj.omega^2*x(3)-2*obj.omega*x(2);
            
            dx = zeros(size(x,1),1);
            dx(1:6:end) = x(2:6:end);
            dx(2:6:end) = a(:,1)+rotTermx;
            dx(3:6:end) = x(4:6:end);
            dx(4:6:end) = a(:,2)+rotTermy;
            dx(5:6:end) = x(6:6:end);
            dx(6:6:end) = a(:,3);
            

        end
        function dx = GoverningEquationsVariational( obj, t, x)
         % variational governing equation in the inertial. 
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
          
            Rot = [cos(obj.omega*t), -sin(obj.omega*t), 0;
                   sin(obj.omega*t),  cos(obj.omega*t), 0;
                   0,                 0,                1];
            
            P = [x(1:6:end),x(3:6:end),x(5:6:end)]*Rot; 
            
            a = obj.gravityModel.acceleration(P);
            gradA = obj.gravityModel.gravityGradient(P);
            
            % correct this ... that not a real tensor
            gradA = (Rot*gradA*Rot');
            a = a*Rot';
            
            dx = zeros(length(x),1);
            dx(1:12:end) = x(2:12:end);
            dx(2:12:end) = a(:,1);
            dx(3:12:end) = x(4:12:end);
            dx(4:12:end) = a(:,2);
            dx(5:12:end) = x(6:12:end);
            dx(6:12:end) = a(:,3);
            dx(7:12:end) = x(8:12:end);
            dx(8:12:end) = gradA(:,1).*x(7 :12:end) + ...
                           gradA(:,2).*x(9 :12:end) + ...
                           gradA(:,3).*x(11:12:end); 
            dx(9:12:end) = x(10:12:end);
            dx(10:12:end) = gradA(:,2).*x(7 :12:end) + ...
                            gradA(:,4).*x(9 :12:end) + ...
                            gradA(:,5).*x(11:12:end); 
            dx(11:12:end) = x(12:12:end);
            dx(12:12:end) = gradA(:,3).*x(7 :12:end) + ...
                            gradA(:,5).*x(9 :12:end) + ...
                            gradA(:,6).*x(11:12:end); 
        end
        
    end
end

