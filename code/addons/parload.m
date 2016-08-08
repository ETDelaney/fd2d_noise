function [y] = parload( filename )

    foo = load(filename);
    whichVariables = fieldnames(foo);

    if numel(whichVariables) == 1
        y = foo.(whichVariables{1});
    else
        error(['Problem with ' filename '!']);
    end

end