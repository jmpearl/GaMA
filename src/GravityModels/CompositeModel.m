classdef CompositeModel
% Composite model 
%--------------------------------------------------------------------------
% acts as a container for multiple models so that SRP or hybrid gravity
% models can be implemented
%--------------------------------------------------------------------------

    properties(GetAccess=public)
        models    % list of models
        frames    % list of model frames (inertial or BFF)
        numModels % number of models
    end    
    methods
        function [obj] = CompositeModel(varargin)
        % initial with set of physics models
        %------------------------------------------------------------------
        % Inputs:
        %   varargin -- allows multiple gravity models to be used in
        %               conjunction
        %------------------------------------------------------------------

            obj.models = varargin;
            obj.numModels = length(varargin);
            
        end
        
        function [obj] = appendModel(obj,model)
        % Provides several methods to distribute mascon throughout interior
        %------------------------------------------------------------------
        % Inputs:
        %   varargin -- allows multiple gravity models to be used in
        %               conjunction
        %------------------------------------------------------------------

            obj.numModels = obj.numModels+1;
            obj.models{obj.numModels} = model;
            
        end
        function [obj] = clearModels(obj)
        % Provides several methods to distribute mascon throughout interior
        %------------------------------------------------------------------
        % Inputs:
        %   varargin -- allows multiple gravity models to be used in
        %               conjunction
        %------------------------------------------------------------------

            obj.numModels = 0;
            obj.models{:} = [];
            
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
               
            potential = zeros(size(p,1),1);
            
            for i = 1:obj.numModels
                
                potential=potential+obj.models{i}.potential(p);
                
            end
            
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
              
            acceleration = zeros(size(p,1),3);

            for i = 1:obj.numModels
                
                acceleration=acceleration+obj.models{i}.acceleration(p);
                
            end
            
        end
        function [laplacian] = laplacian(obj,p)
        % Laplacian of gravitational field. Always zero for point mass 
        %------------------------------------------------------------------
        % Inputs:
        %   p --------- [Mx3] array of M sample coordinates
        %------------------------------------------------------------------
        % Outputs:
        %   laplacian - [Mx1] laplacian at sample sites
        %------------------------------------------------------------------
          
            laplacian = zeros(size(p,1),1);
            
            for i = 1:obj.numModels
                
                laplacian=laplacian+obj.models{i}.laplacian(p);
                
            end
        end
        function [gravGradient] = gravityGradient(obj,p)
        % Symmetric gravitational gradient tensor stored as a 1x6
        %------------------------------------------------------------------
        % Inputs:
        %   p ------------ [Mx3] array of M sample coordinates
        %------------------------------------------------------------------
        % Outputs:
        %   gravGradient - [Mx6] acceleration gradient at sample sites
        %------------------------------------------------------------------
         
            gravGradient = zeros(size(p,1),6);
            for i = 1:obj.numModels
                
                gravGradient=gravGradient+obj.models{i}.gravGradient(p);
                
            end

        end
    end
end
