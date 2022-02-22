function [isContained] = contains(string,substring)
    isContained = ~isempty(strfind(string,substring));
end