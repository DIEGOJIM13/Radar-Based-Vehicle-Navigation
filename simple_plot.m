%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Viasat Radar Based Behicle Location and Navigation System
%University of Arizona ENG498 Team 16060

%Simple_Plot Software
%Comment:
    %This file plots the time domain inputted sets of IQ data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;
fileNum = sprintf('%d',3);
path1 = 'data/R1_5-4_set3/';
path2 = 'data/R2_5-4_set3/';
filename1 = strcat(path1, fileNum, '.txt');
filename2 = strcat(path2, fileNum, '.txt');
data1 = load(filename1);
data2 = load(filename2);
length1 = size(data1);
length2 = size(data2);

t1 = data1(1:length1,1);
I1 = data1(1:length1,2);
Q1 = data1(1:length1,3);
t2 = data2(1:length2,1);
I2 = data2(1:length2,2);
Q2 = data2(1:length2,3);

plot(t1,I1,'b-','Markersize',15);
hold on
plot(t1,Q1,'g-','Markersize',15);
title('Radar1 (Passenger) Sampled IQ Data');
xlabel('n');
ylabel('Voltage (V)');
legend('I','Q');

figure;
plot(t2,I2,'b-','Markersize',15);
hold on
plot(t2,Q2,'g-','Markersize',15);
title('Radar2 (Driver) Sampled IQ Data');
xlabel('n');
ylabel('Voltage (V)');
legend('I','Q');

figure;
plot3(t1,I1,Q1,'b.');
title('Radar1 (Passenger)');
xlabel('n');
ylabel('Voltage I');
zlabel('Voltage Q');

figure;
plot3(t2,I2,Q2,'b.');
title('Radar2 (Driver)');
xlabel('n');
ylabel('Voltage I');
zlabel('Voltage Q');

%meanI1 = mean(I1)
%meanQ1 = mean(Q1)
%meanI2 = mean(I2)
%meanQ2 = mean(Q2)

if(length1 ~= length2)
    display('lengths do not equal')
end