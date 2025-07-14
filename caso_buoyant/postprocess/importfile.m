function spectra2 = importdata(filename, startRow, endRow)

%Initialize variables.
if nargin<=2
    startRow = 3;
    endRow = inf;
end

%Format for each line of text:
formatSpec = '%24f%25f%25f%25f%25f%25f%25f%25f%f%[^\n\r]';

%Open the text file.
fileID = fopen(filename,'r');

%Read columns of data according to the format.
% This call is based on the structure of the file used to generate this code. 
textscan(fileID, '%[^\n\r]', startRow(1)-1, 'WhiteSpace', '', 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    textscan(fileID, '%[^\n\r]', startRow(block)-1, 'WhiteSpace', '', 'ReturnOnError', false);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%Close the text file.
fclose(fileID);

%Create output variable
spectra2 = [dataArray{1:end-1}];
