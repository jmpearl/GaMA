classdef  progressbar < handle
% simple class to monitor the progress of a calculation
%--------------------------------------------------------------------------
% OOP implementation of textprogressbar by Paul Proteus
% https://www.mathworks.com/matlabcentral/fileexchange/28067-text-progress-bar
%--------------------------------------------------------------------------
    properties (GetAccess=public)
        strCR;   
        strPercentageLength = 10;  
        strDotsMaximum      = 10;  
    end

    methods (Access=public)
        function obj = progressbar(str)
            if nargin==1
                obj.initialize(str);
            end
        end
        function initialize(obj,str)
            fprintf('%s',str);
            obj.strCR=-1;
        end
        function update(obj,num)
            
            num = floor(num);
            percentageOut = [num2str(num) '%%'];
            percentageOut = [percentageOut repmat(' ',1,obj.strPercentageLength-length(percentageOut)-1)];
            nDots = floor(num/100*obj.strDotsMaximum);
            dotOut = ['[' repmat('.',1,nDots) repmat(' ',1,obj.strDotsMaximum-nDots) ']'];
            strOut = [percentageOut dotOut];

            if obj.strCR == -1
                fprintf(strOut);
            else
                fprintf([obj.strCR strOut]);
            end

            % Update carriage return
            obj.strCR = repmat('\b',1,length(strOut)-1);
        end
        function close(obj)
            obj.update(100);
            fprintf(' done.\n');
        end
    end
end
