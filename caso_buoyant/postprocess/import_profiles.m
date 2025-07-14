function profiles = import_profiles(filename)
    % Open the file
    fid = fopen(filename, 'r');
    if fid == -1
        error('Cannot open file: %s', filename);
    end
    
    % Read all lines
    lines = {};
    tline = fgetl(fid);
    while ischar(tline)
        lines{end+1} = tline; %#ok<AGROW>
        tline = fgetl(fid);
    end
    fclose(fid);
    
    % First pass: determine number of columns
    ncols = 0;
    for i = 1:length(lines)
        line = strtrim(lines{i});
        if isempty(line) || startsWith(line, '#')
            continue
        end
        nums = sscanf(line, '%f');
        if ~isempty(nums)
            ncols = length(nums);
            break
        end
    end
    
    if ncols == 0
        error('Could not determine the number of columns: no valid data lines found.');
    end

    % Now parse each line
    data = [];
    for i = 1:length(lines)
        line = strtrim(lines{i});
        if isempty(line) || startsWith(line, '#')
            continue  % Skip comment lines
        else
            nums = sscanf(line, '%f')';
            % Keep only lines with exactly ncols numbers
            if length(nums) == ncols
                data = [data; nums];
            end
        end
    end
    
    profiles = data;
end
