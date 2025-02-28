function   [r,v]=readDat(datFileName,headerLength)
fid = fopen(datFileName,'r');           
charData = [''];
lineCounter = 1;
tline = fgetl(fid);
while ischar(tline)
    str = strtrim(tline);
    str = strrep(str,',',' ');
    %disp(str)
    if ~isempty(tline) && lineCounter >= headerLength
        charData = [charData,str];
    end
    lineCounter = lineCounter+1;
    tline = fgetl(fid);
end
fclose(fid);
data = cell2mat(textscan(charData,'%f%f%f%f%f%f'));
r = data(:,1:3);
v = data(:,4:6);
end