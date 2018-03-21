%% 2018 ViaSat Radar Navigation Post Processing Algorithm Team 1718 
% This algorithm takes in raw radar data Is and Qs and converts them into
%doplar frequencies, then to velocities and finally to positions. In order
%to do this accurately this algorithm implements serveral filters that
%denoise the radar data.

%% Clear Stuff
clear all;
close all;

%%  Define Variables
fs = 18000; %Hz hardware sampling rate
fc_d = 24160000000; %Hz carrier freq driver side 
fc_p = 24150000000; %Hz carrier freq passanger side
c_vac = 299792458; %m/s speed of light in vacuum
c_air = c_vac/1.0003; %speed of light in air
t = 1/fs; %time per sample

%% Read in Radar Data
filename1 = strcat('RealTest.txt'); %open driver side txt file
filename2 = strcat('RealTest.txt'); %open passanger side txt file

data = load(filename1); %load driver side
length = size(data); %find length
I1 = data(1:length,2)'; %get In phase data- driver
Q1 = data(1:length,3)'; %get Quadrature data-driver
    
data2 = load(filename2); %load passanger data 
length = size(data); %find length
I2 = data2(1:length,2)'; %get In phase data- passanger
Q2 = data2(1:length,3)';%get Quadrature data- passanger

%% Recenter I and Q data and low pass
 %{   
I1_t = I1 - mean(I1);
Q1_t = Q1 - mean(Q1);
I2_t = I2 - mean(I2);
Q2_t = Q2 - mean(Q2);

fcut = 5000;
[b,a] = butter(6,fcut/(fs/2));
I1 = filter(b,a,I1_t);
Q1 = filter(b,a,Q1_t);
I2 = filter(b,a,I2_t);
Q2 = filter(b,a,Q2_t);

%}

%% convert to complex
dataCplx1 = complex(I1,Q1);
absArray1 = abs(dataCplx1);
phaseArray1 = angle(dataCplx1); %find phase - driver 

dataCplx2 = complex(I2,Q2);
absArray2 = abs(dataCplx2);
phaseArray2 = angle(dataCplx2);%find phase - passanger

%% convert to doppler frequency

fd = diff(phaseArray1)/(2*pi); %driver doplar frequency
fp = diff(phaseArray2)/(2*pi); %passanger doplar frequency

%% convert to velocity 

vd = (fd*c_air)/(2*fc_d);%find velocity driver
vp = (fp*c_air)/(2*fc_p);%find velocity passanger 

%% calculate distance 

dd = vd*t;% distance vd
dp = vp*t;

dtd = cumsum(dd);
dtp = cumsum(dp);
%disp(dt)


%% plot path

%size = s;
%xlims = [-size size]; ylims = [-size size];
%line([xlims nan 0 0],[0 0 nan ylims],'LineWidth',0.5, 'Color',[.2 .2 .2])
%axis square, grid on
%set(gca, 'XLim',xlims, 'YLim',ylims)
%title('Vehicle Path')
%xlabel('x distance (m)'), ylabel('y distance (m)')

%hold on

%plot(1,dtv, 'LineStyle','none', 'Marker','o', 'Color','m');


