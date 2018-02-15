%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Viasat Radar Based Behicle Location and Navigation System
%University of Arizona ENG498 Team 16060

%Bin2txt Software
%Comment:
    %This File reads in binary files and converts them into newly created
    %text files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
tic;
for fileNum = 8:8%change these number to loop through all binary files
    fileNum_char = sprintf('%d', fileNum);
    filename = strcat('data/R2_5-4_set3/',fileNum_char);%change input file folder

    fileID = fopen(strcat(filename,'.bin'),'r');
    A = fread(fileID,'uint16');
    fclose(fileID);
    
    n = 1;
    i = 1;
    j = 1;
    %convert 12bit ADC values to voltages
    while i<=size(A,1)
        data(j,1) = n;
        data(j+1,1) = A(i)*3.3/4095;
        data(j+2,1) = A(i+1)*3.3/4095;
        n = n + 1;
        i = i + 2;
        j = j + 3;
    end
    
    %write to txt file
    fileID = fopen(strcat(filename,'.txt'),'w');
    fprintf(fileID, '%07d\t%.3f\t%.3f\n', data);
    fclose(fileID);
end
toc;