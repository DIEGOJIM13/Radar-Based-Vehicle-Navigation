%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Viasat Radar Based Behicle Location and Navigation System
%University of Arizona ENG498 Team 16060

%Bin2txt Software
%Comment:
    %This File reads in binary files and converts them into newly created
    %text files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [left, right] = bin2txt(fileName1, fileName2)
    if (exist('fileName1','var') && exist('fileName2','var'))
        leftF = fileName1;
        rightF = fileName2;
    elseif (exist('fileName1','var'))
        leftF = fileName1;
        rightF = fileName1;
    else
        left = []; right = [];
        [filenameLeft, pathnameLeft] = uigetfile('*.bin', 'Select the left .bin file');
        if isequal(filenameLeft,0)
           return;
        else
           disp(['Left File  ', fullfile(pathnameLeft, filenameLeft)]);
        end
        [filenameRight, pathnameRight] = uigetfile('*.bin', 'Select the right .bin file');
        if isequal(filenameRight,0)
           return;
        else
           disp(['Right File ', fullfile(pathnameRight, filenameRight)]);
        end
        leftF = fullfile(pathnameLeft, filenameLeft);
    	rightF = fullfile(pathnameRight, filenameRight);
    end
    
    leftFileID = fopen(leftF, 'r');
    leftRaw = fread(leftFileID, 'uint16');
    fclose(leftFileID);
    
    
    rightFileID = fopen(rightF, 'r');
    rightRaw = fread(rightFileID, 'uint16');
    fclose(rightFileID);

    left = zeros(fix(size(leftRaw, 1) / 2), 3);
    right = zeros(fix(size(rightRaw, 1) / 2), 3);
    
    i = 1; n = 1;
    %convert 12bit ADC values to voltages
    while i <= fix(size(leftRaw, 1) / 2)
        left(i,1) = i;
        left(i,2) = leftRaw(n)*3.3/4095;
        left(i,3) = leftRaw(n+1)*3.3/4095;
        i = i + 1; n = n + 2;
    end

    i = 1; n = 1;
    while i <= fix(size(rightRaw, 1) / 2)
        right(i,1) = i;
        right(i,2) = rightRaw(n)*3.3/4095;
        right(i,3) = rightRaw(n+1)*3.3/4095;
        i = i + 1; n = n + 2;
    end
    

