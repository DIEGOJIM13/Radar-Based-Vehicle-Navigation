%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Viasat Radar Based Behicle Location and Navigation System
%University of Arizona ENG498 Team 16060

%Data Analysis Software
%Comment:
    %This file preforms phase difference and FFT analysis on the input IQ
    %data and returns a velocity, location, and heading
    %need to add heading
    %need to add better numerical integration
    %need to add kalman filter/ data filtering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

%variables
fs = 18000; %Hz hardware sampling rate
fc_driver = 24160000000; %Hz carrier freq
fc_passenger = 24150000000; %Hz 
c_vac = 299792458; %m/s speed of light in vacuum
c_air = c_vac/1.0003; %speed of light in air
%[t, I1, Q1] = constSpeedIQgen(9.872129045*sin(pi/4),fs,fc,c_air,30); %generate IQ data
%[t, I1, Q1] = constSpeedIQgen(10.17926195*sin(pi/4),fs,fc,c_air,30); %generate IQ data
%[t, I1, Q1] = constSpeedIQgen(35*sin(pi/4),fs,fc,c_air,5);
%[t, I1, Q1] = IQaccel(40*sin(pi/4),fs,fc,c_air,5); %generate IQ data
fileNum = 3;
fileNum_char = sprintf('%d', fileNum);
filename1 = strcat('tests/R1_5-4_set3/',fileNum_char,'.txt');
filename2 = strcat('tests/R2_5-4_set3/',fileNum_char,'.txt');
[t1, I1, Q1] = readIQ(filename1);
[t2, I2, Q2] = readIQ(filename2);

I1 = I1 - mean(I1);
Q1 = Q1 - mean(Q1);
I2 = I2 - mean(I2);
Q2 = Q2 - mean(Q2);

%Low pass butterworth filter
fcut = 5000;
[b,a] = butter(6,fcut/(fs/2));
I1 = filter(b,a,I1);
Q1 = filter(b,a,Q1);
I2 = filter(b,a,I2);
Q2 = filter(b,a,Q2);

figure('Name','Radar1 (Passenger) IQplot');
plot3(t1, I1, Q1, '.'); %plot IQ data

figure('Name','Radar2 (Driver) IQplot');
plot3(t2, I2, Q2, '.'); %plot IQ data

%convert IQ to complex form
dataCplx1 = complex(I1,Q1);
absArray1 = abs(dataCplx1);
phaseArray1 = angle(dataCplx1);

%figure;
%spectrogram(dataCplx1,256,250,256,fs,'yaxis')

dataCplx2 = complex(I2,Q2);
absArray2 = abs(dataCplx2);
phaseArray2 = angle(dataCplx2);

%determine phase difference to determine forward or backward movement
tlength = min(length(t1),length(t2));
t = (1/fs)*t1(1,1:tlength);
phaseDiff1 = zeros(1,length(t));
phaseDiff2 = zeros(1,length(t));
for n = 2:length(t)
    phaseDiff1(n) = mod(phaseArray1(n) - phaseArray1(n-1) + 3*pi, 2*pi) - pi;
    phaseDiff2(n) = mod(phaseArray2(n) - phaseArray2(n-1) + 3*pi, 2*pi) - pi;
end
phaseDiff1(1) = phaseDiff1(2); %assume first sample is same as second (for now)
phaseDiff2(1) = phaseDiff2(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Estimations with simple phase difference estimate
%estimate doppler freqs and velocities
% fdopp1_pd = zeros(1,length(t));
% fdopp2_pd = zeros(1,length(t));
% vx_pd = zeros(1,length(t));
% vy_pd = zeros(1,length(t));
% distDiffx_pd = zeros(1,length(t));
% distDiffy_pd = zeros(1,length(t));
% distTotx_pd = zeros(1,length(t));
% distToty_pd = zeros(1,length(t));
% 
% distDiffx_pd_mean = zeros(1,length(t));
% distDiffy_pd_mean = zeros(1,length(t));
% distTotx_pd_mean = zeros(1,length(t));
% distToty_pd_mean = zeros(1,length(t));
% 
% n = 2;
% [fdopp1_pd(n), fdopp2_pd(n), vx_pd(n), vy_pd(n)] = estVelocity_PD(2, phaseDiff1, phaseDiff2, fs, c_air, fc_driver, fc_passenger, vx_pd, vy_pd);
% vx_pd(2) = 2*vx_pd(2);
% vy_pd(2) = 2*vy_pd(2);
% fdopp1_pd(1) = fdopp1_pd(2);
% fdopp2_pd(1) = fdopp2_pd(2);
% for n = 3:length(t)
%     [fdopp1_pd(n), fdopp2_pd(n), vx_pd(n), vy_pd(n)] = estVelocity_PD(n, phaseDiff1, phaseDiff2, fs, c_air, fc_driver, fc_passenger, vx_pd, vy_pd);
% end
% k = 18000;
% vx_pd_mean = movmean(vx_pd,k);
% vy_pd_mean = movmean(vy_pd,k);
% figure;
% plot(t,vx_pd_mean,'b.');
% hold on
% plot(t,vy_pd_mean,'r.');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%estimations with FFT
N = 1024; %N-point FFT
Z = 16/4; %Zeropadding multiplier
fdopp1_fft = zeros(1,length(t));
fdopp2_fft = zeros(1,length(t));
vx_fft = zeros(1,length(t));
vy_fft = zeros(1,length(t));
distDiffx_fft = zeros(1,length(t));
distDiffy_fft = zeros(1,length(t));
distTotx_fft = zeros(1,length(t));
distToty_fft = zeros(1,length(t));

distDiffx_fft_mean = zeros(1,length(t));
distDiffy_fft_mean = zeros(1,length(t));
distTotx_fft_mean = zeros(1,length(t));
distToty_fft_mean = zeros(1,length(t));

%FFT for first N samples
n = N;
dft1 = fft(dataCplx1((n-(N-1)):n),Z*N);
[peak1,freqInd1] = max(abs(dft1(1,1:length(dft1)/2)));
dft2 = fft(dataCplx2((n-(N-1)):n),Z*N);
[peak2,freqInd2] = max(abs(dft2(1,1:length(dft2)/2)));
freqArr = 0:fs/(Z*N):fs/2;

%Filter out no movement/low magnitude noise
if(peak1 > 5)
    fdopp1_fft((n-(N-1)):n) = freqArr(freqInd1)*cos(.123*pi);
else
    fdopp1_fft((n-(N-1)):n) = 0;
end
if(peak2 > 5)
    fdopp2_fft((n-(N-1)):n) = freqArr(freqInd2)*cos(.123*pi);
else
    fdopp2_fft((n-(N-1)):n) = 0;
end

%Convert Doppler frequencies into velocities
vx_fft((n-(N-1)):n) = ((fdopp2_fft((n-(N-1)):n)/fc_driver+fdopp1_fft((n-(N-1)):n)/fc_passenger)*c_air/(2)*sin(pi/4));
vy_fft((n-(N-1)):n) = ((fdopp2_fft((n-(N-1)):n)/fc_driver-fdopp1_fft((n-(N-1)):n)/fc_passenger)*c_air/(2)*cos(pi/4));

%Loop through using 1 sample overlapping FFT
for n = (N+1):length(t)
    dft1 = fft(dataCplx1((n-(N-1)):n),Z*N);
    [peak1,freqInd1] = max(abs(dft1(1,1:length(dft1)/2)));
    dft2 = fft(dataCplx2((n-(N-1)):n),Z*N);
    [peak2,freqInd2] = max(abs(dft2(1,1:length(dft2)/2)));
    
    %Filter out no movement/low magnitude noise
    if(peak1 > 5)
        fdopp1_fft((n-(N-1)):(n-1)) = (fdopp1_fft((n-(N-1)):(n-1)) + freqArr(freqInd1)*cos(.123*pi))/2;
        fdopp1_fft(n) = freqArr(freqInd1)*cos(.123*pi);
    else
        fdopp1_fft((n-(N-1)):(n-1)) = (fdopp1_fft((n-(N-1)):(n-1)) + 0)/2;
        fdopp1_fft(n) = 0;
    end
    if(peak2 > 5)
        fdopp2_fft((n-(N-1)):(n-1)) = (fdopp2_fft((n-(N-1)):(n-1)) + freqArr(freqInd2)*cos(.123*pi))/2;
        fdopp2_fft(n) = freqArr(freqInd2)*cos(.123*pi);
    else
        fdopp2_fft((n-(N-1)):(n-1)) = (fdopp2_fft((n-(N-1)):(n-1)) + 0)/2;
        fdopp2_fft(n) = 0;
    end
    
    %Convert Doppler frequencies into velocities
    vx_fft((n-(N-1)):(n-1)) = ((fdopp2_fft((n-(N-1)):(n-1))/fc_driver+fdopp1_fft((n-(N-1)):(n-1))/fc_passenger)*c_air/(2)*sin(pi/4));
    vx_fft(n) = ((fdopp2_fft(n)/fc_driver+fdopp1_fft(n)/fc_passenger)*c_air/(2)*sin(pi/4));
    vy_fft((n-(N-1)):(n-1)) = ((fdopp2_fft((n-(N-1)):(n-1))/fc_driver-fdopp1_fft((n-(N-1)):(n-1))/fc_passenger)*c_air/(2)*cos(pi/4));
    vy_fft(n) = ((fdopp2_fft(n)/fc_driver-fdopp1_fft(n)/fc_passenger)*c_air/(2)*cos(pi/4));
end

k = 18000;
vx_fft_mean = movmean(vx_fft,k);
vy_fft_mean = movmean(vy_fft,k);
figure;
plot(t,vx_fft_mean,'b.');
hold on
plot(t,vy_fft_mean,'r.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%distance estimation
% n = 5;
% distDiffx_pd(n) = (4/90)*(7*vx_pd(n-4)+32*vx_pd(n-3)+12*vx_pd(n-2)+32*vx_pd(n-1)+7*vx_pd(n))*(1/fs);
% distDiffy_pd(n) = (4/90)*(7*vy_pd(n-4)+32*vy_pd(n-3)+12*vy_pd(n-2)+32*vy_pd(n-1)+7*vy_pd(n))*(1/fs);
% distTotx_pd(n) = distDiffx_pd(n);
% distToty_pd(n) = distDiffy_pd(n);
% 
% distDiffx_fft(n) = (4/90)*(7*vx_fft(n-4)+32*vx_fft(n-3)+12*vx_fft(n-2)+32*vx_fft(n-1)+7*vx_fft(n))*(1/fs);
% distDiffy_fft(n) = (4/90)*(7*vy_fft(n-4)+32*vy_fft(n-3)+12*vy_fft(n-2)+32*vy_fft(n-1)+7*vy_fft(n))*(1/fs);
% distTotx_fft(n) = distDiffx_fft(n);
% distToty_fft(n) = distDiffy_fft(n);
% for n = 10:5:length(t)-mod(length(t),5)

%Straight Line Approximations
for n = 2:length(t)

%     distDiffx_pd(n) = .5*(vx_pd(n) + vx_pd(n-1))*(1/fs);
%     distDiffy_pd(n) = .5*(vy_pd(n) + vy_pd(n-1))*(1/fs);
%     distTotx_pd(n) = distTotx_pd(n-1) + distDiffx_pd(n);
%     distToty_pd(n) = distToty_pd(n-1) + distDiffy_pd(n);
    
    %FFT straight line distance approximations
    distDiffx_fft(n) = .5*(vx_fft(n) + vx_fft(n-1))*(1/fs);
    distDiffy_fft(n) = .5*(vy_fft(n) + vy_fft(n-1))*(1/fs);
    distTotx_fft(n) = distTotx_fft(n-1) + distDiffx_fft(n);
    distToty_fft(n) = distToty_fft(n-1) + distDiffy_fft(n);
    
%     distDiffx_pd_mean(n) = .5*(vx_pd_mean(n) + vx_pd_mean(n-1))*(1/fs);
%     distDiffy_pd_mean(n) = .5*(vy_pd_mean(n) + vy_pd_mean(n-1))*(1/fs);
%     distTotx_pd_mean(n) = distTotx_pd_mean(n-1) + distDiffx_pd_mean(n);
%     distToty_pd_mean(n) = distToty_pd_mean(n-1) + distDiffy_pd_mean(n);
    
    %Moving 2 point average of distance
    distDiffx_fft_mean(n) = .5*(vx_fft_mean(n) + vx_fft_mean(n-1))*(1/fs);
    distDiffy_fft_mean(n) = .5*(vy_fft_mean(n) + vy_fft_mean(n-1))*(1/fs);
    distTotx_fft_mean(n) = distTotx_fft_mean(n-1) + distDiffx_fft_mean(n);
    distToty_fft_mean(n) = distToty_fft_mean(n-1) + distDiffy_fft_mean(n);

    %Boole's rule approximation (Newton-Cotes Formula degree 4)
%     distDiffx_pd(n) = (4/90)*(7*vx_pd(n-4)+32*vx_pd(n-3)+12*vx_pd(n-2)+32*vx_pd(n-1)+7*vx_pd(n))*(1/fs);
%     distDiffy_pd(n) = (4/90)*(7*vy_pd(n-4)+32*vy_pd(n-3)+12*vy_pd(n-2)+32*vy_pd(n-1)+7*vy_pd(n))*(1/fs);
%     distTotx_pd(n) = distTotx_pd(n-5) + distDiffx_pd(n);
%     distToty_pd(n) = distToty_pd(n-5) + distDiffy_pd(n);
%     
%     distDiffx_fft(n) = (4/90)*(7*vx_fft(n-4)+32*vx_fft(n-3)+12*vx_fft(n-2)+32*vx_fft(n-1)+7*vx_fft(n))*(1/fs);
%     distDiffy_fft(n) = (4/90)*(7*vy_fft(n-4)+32*vy_fft(n-3)+12*vy_fft(n-2)+32*vy_fft(n-1)+7*vy_fft(n))*(1/fs);
%     distTotx_fft(n) = distTotx_fft(n-5) + distDiffx_fft(n);
%     distToty_fft(n) = distToty_fft(n-5) + distDiffy_fft(n);
end

% figure('Name','Velocity_PD');
% plot(t,vx_pd,'b.');
% hold on
% plot(t,vy_pd,'r.');
% title('Velocity Estimation Using Phase Difference');
% xlabel('time (s)');
% ylabel('velocity (m/s)');

% figure('Name','Location_PD');
% plot(distTotx_pd,distToty_pd,'.');
% title('Location Estimation Using Phase Difference');
% xlabel('x distance (m)');
% ylabel('y distance (m)');

%Plots
%FFT Velocity Plot
figure('Name','Velocity_FFT');
plot(t,vx_fft,'b.');
hold on
plot(t,vy_fft,'r.');
title('Velocity Estimation Using FFT');
xlabel('time (s)');
ylabel('velocity (m/s)');

%FFT Location Plot
figure('Name','Location_FFT');
plot(distTotx_fft,distToty_fft,'.');
title('Location Estimation Using FFT');
xlabel('x distance (m)');
ylabel('y distance (m)');

%FFT Location Plot with Moving Average
figure('Name','Location_FFT');
plot(distTotx_fft_mean,distToty_fft_mean,'.');
title('Location Estimation Using FFT');
xlabel('x distance (m)');
ylabel('y distance (m)');

toc;
%FFT or Radar Return
figure;
plot(freqArr, abs(dft1(1:(Z*N)/2+1)),'-');
title('FFT of Radar Doppler Return');
xlabel('frequency (Hz)');
ylabel('Magnitude');