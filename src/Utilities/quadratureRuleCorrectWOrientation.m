directory = 'quadratureSchemes/prism/';
files = dir(directory);
for i = 1:length(files)
    if ~files(i).isdir
        name = files(i).name;
        load([directory,name]);
        quadratureScheme.weights = quadratureScheme.weights';

        save([directory,name],'quadratureScheme')
    end
end

