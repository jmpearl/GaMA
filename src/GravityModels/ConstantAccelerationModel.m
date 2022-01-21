classdef ConstantAccelerationModel
% applies constant acceleration in user defined direction 
%--------------------------------------------------------------------------
% for third body gravity fields, sun, jupiter etc...
%--------------------------------------------------------------------------

    properties(GetAccess=public)
        frame; % integrator needs to know if bff or inertial x,y,z
        a;     % acceleration
    end    
    methods
        function [obj] = ConstantAccelerationModel(a)
        % initial with set of physics models
        %------------------------------------------------------------------
        % Inputs:
        %   a -- constant acceleration
        %------------------------------------------------------------------
            obj.frame = 'inertial';
            obj.a = a;
            
        end
        function [obj] = addAcceleration(obj,a)
           obj.a = obj.a + a; 
        end
        function [potential] = potential(obj,p)
        % Gravitational potential using the negative convention
        %------------------------------------------------------------------
        % Inputs:
        %   p --------- Mx3 array of M sample coordinates
        %------------------------------------------------------------------
        % Outputs:
        %   potential - Mx1 gravitational potential at sample sites
        %------------------------------------------------------------------
               
            potential = p*obj.a';
            
        end
        function [acceleration] = acceleration(obj,p)
        % Gravitational acceleration 
        %------------------------------------------------------------------
        % Inputs:
        %   p ------------ [Mx3] array of M sample coordinates
        %------------------------------------------------------------------
        % Outputs:
        %   acceleration - [Mx3] gravitational acceleration at sample sites
        %------------------------------------------------------------------
              
            acceleration = zeros(size(p,1),3)+obj.a;
            
        end
    end
    methods(Static,Access=public)
        function [laplacian] = laplacian(p)
        % Laplacian of gravitational field. Always zero for point mass 
        %------------------------------------------------------------------
        % Inputs:
        %   p --------- [Mx3] array of M sample coordinates
        %------------------------------------------------------------------
        % Outputs:
        %   laplacian - [Mx1] laplacian at sample sites
        %------------------------------------------------------------------
          
            laplacian = zeros(size(p,1),1);
            
        end
        function [gravGradient] = gravityGradient(p)
        % Symmetric gravitational gradient tensor stored as a 1x6
        %------------------------------------------------------------------
        % Inputs:
        %   p ------------ [Mx3] array of M sample coordinates
        %------------------------------------------------------------------
        % Outputs:
        %   gravGradient - [Mx6] acceleration gradient at sample sites
        %------------------------------------------------------------------
         
            gravGradient = zeros(size(p,1),6);

        end
    end
end
