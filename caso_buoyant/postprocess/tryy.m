function profiles = import_profiles(filename, startRow)
%IMPORT_PROFILES Import data with possible non-numeric entries.

if nargin < 2
    startRow = 3;
end

% Define number of columns
ncols = 34;

fileID = fopen(filename, 'r');
if fileID == -1
    error('Could not open file.');
end

% Read everything as strings
dataArray = textscan(fileID, repmat('%s',1,ncols), ...
    'HeaderLines', startRow-1, ...
    'Delimiter', '', ...
    'MultipleDelimsAsOne', true, ...
    'ReturnOnError', false);

fclose(fileID);

% Now process manually:
nrows = length(dataArray{1});
profiles = cell(nrows, ncols); % cell array: strings or numbers

for i = 1:ncols
    for j = 1:nrows
        str = dataArray{i}{j};
        num = str2double(str);
        if isnan(num) && ~strcmpi(str, 'NaN')
            % Keep original string if it is non-numeric and not literal NaN
            profiles{j,i} = str;
        else
            % Otherwise store the numeric value
            profiles{j,i} = num;
        end
    end
end
end
