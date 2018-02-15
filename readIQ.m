%readIQ function reads in input IQ data
function [ t, I, Q ] = readIQ( filename )
    data = load(filename);
    length = size(data);
    t = data(1:length,1)';
    I = data(1:length,2)';
    Q = data(1:length,3)';
end

