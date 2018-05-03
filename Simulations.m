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
c_vac = 299792458; %m/s speed of light in vacuum
t = 1/fs; %time per sample
d = 0.5; %dist between radar modules 
rng(1);

%% Get Noise Data from Previous Data 
filename1 = strcat('Book1.txt'); %open driver side txt file
filename2 = strcat('Book1.txt'); %open passanger side txt file

data =load(filename1); %load driver side
%data = a;
Id = data(1:size(data),2)'; %get In phase data- driver
Qd = data(1:size(data),3)'; %get Quadrature data-driver

dataCplxd = complex(Id,Qd);
phaseArrayd = angle(dataCplxd);

fd_dop = fs*(diff(phaseArrayd)/(pi));
v_d = (fd_dop*c_vac)/(2*fc_d);

Amp = v_d(1:1.5*fs); 

Mag = sqrt(Id.^2 + Qd.^2);
avg = mean(Mag);
stand = std(Mag);
snr = 20*log(avg/stand);


%% Create Radar Data

sec = 5;
noise = 60; % snr dB
vel = 30;
n = 198000;

i = int64(n/4);
lamdah1 = c_vac/fc_d;

% v_start = (1:(18000*5))/18000;
% v_end = (1:(18000*6))/18000;
% v_end = v_end + max(v_start) / 0.1;
% v_end = 0.1 * v_end;
% v1 = [v_start, v_end];
% v2 = [1.01 * v_start, v_end];

v1 = [2.0*ones(18000*3, 1); 3.0*ones(18000*2.4, 1); 3.0*ones(18000*0.6, 1); 4.0*ones(18000, 1); 5.0*ones(18000, 1);  6.0*ones(18000, 1);  7.0*ones(18000, 1);  7.5*ones(18000, 1)];
v2 = [2.1*ones(18000*3, 1); 3.2*ones(18000*2.4, 1); 3.0*ones(18000*0.6, 1); 4.0*ones(18000, 1); 5.0*ones(18000, 1);  6.0*ones(18000, 1);  7.0*ones(18000, 1);  7.5*ones(18000, 1)];
  
t1 = 1:length(v1);
t1 = (t1*(1/fs))';

fd = speed2dop(v1, lamdah1);
fp = speed2dop(v2, lamdah1);

I1 = cos(2*pi*t1.*fd);
Q1 = sin(2*pi*t1.*fd);

I2 = cos(2*pi*t1.*fp);
Q2 = sin(2*pi*t1.*fp);

I1_n = awgn(I1,noise);
Q1_n = awgn(Q1,noise);

I2_n = awgn(I2,noise);
Q2_n = awgn(Q2,noise);


%% convert to complex

%without noise
dataCplx1 = complex(I1,Q1);
absArray1 = abs(dataCplx1);
phaseArray1 = angle(dataCplx1); %find phase - driver 
phaseArray1 = unwrap(phaseArray1);


dataCplx2 = complex(I2,Q2);
absArray2 = abs(dataCplx2);
phaseArray2 = angle(dataCplx2);%find phase - passanger
phaseArray2 = unwrap(phaseArray2);

%without noise 
dataCplx1_n = complex(I1_n,Q1_n);
absArray1_n = abs(dataCplx1_n);
phaseArray1_n = angle(dataCplx1_n); %find phase - driver 
phaseArray1_n = unwrap(phaseArray1_n);


dataCplx2_n = complex(I2_n,Q2_n);
absArray2_n = abs(dataCplx2_n);
phaseArray2_n = angle(dataCplx2_n);%find phase - passanger
phaseArray2_n = unwrap(phaseArray2_n);



%% convert to doppler frequency

%doppler without noise
fd1 = fs*(diff(phaseArray1)/(pi)); %driver doplar frequency
fp1 = fs*(diff(phaseArray2)/(pi)); %passanger doplar frequency

%doppler with noise
fd1_n = fs*(diff(phaseArray1_n)/(pi)); %driver doplar frequency
fp1_n = fs*(diff(phaseArray2_n)/(pi)); %passanger doplar frequency


%% convert to velocity 
% radial velocity scaling by trig azimuth and elevation angle 

%velocity without noise
vd = (fd1*c_vac)/(2*fc_d); %find velocity driver
%     *cos(sqrt(2)/2)* sin(theta - 90) * -1; % scale for y and z and direction if car moving in x
vp = (fp1*c_vac)/(2*fc_d); %find velocity driver
%     *cos(sqrt(2)/2)* sin(theta - 90) * -1; % scale for y and z and direction if car moving in x

%velocity with noise
%velocity without noise
vd_n = (fd1_n*c_vac)/(2*fc_d); %find velocity driver
%     *cos(sqrt(2)/2)* sin(theta - 90) * -1; % scale for y and z and direction if car moving in x
vp_n = (fp1_n*c_vac)/(2*fc_d); %find velocity driver
%     *cos(sqrt(2)/2)* sin(theta - 90) * -1; % scale for y and z and direction if car moving in x

%% noise filtering

%vd3 = movmean(vd,200);
%vp= movmean(vp,200);

a = 1;
b = [1/4 1/4 1/4 1/4];

%wiener_d = imgaussfilt(vd_n);
%wiener_p = imgaussfilt(vp_n);


wiener_d = wiener2(vd_n, [1, 1000])+mean(vd_n);
wiener_p = wiener2(vp_n, [1, 1000])+mean(vp_n);
 %wiener_d = vd;
 %wiener_p = vp;

%figure(1);
%plot(vd);
%figure(2);
%plot(wiener_d(100:end-100));
%figure(3);
%plot(vd3);



%% calculate angle

%angle without noise
omega =  - (vd - vp) ./ d;

rel_angle = omega * t;

angle = cumsum(rel_angle);

%angle with noise 
omega_n =  - (vd_n - vp_n) ./ d;

rel_angle_n = omega_n * t;

angle_n = cumsum(rel_angle_n);


%angle with noise + filtering 
omega_nf =  - (wiener_d - wiener_p) ./ d;

rel_angle_nf = omega_nf * t;

angle_nf = cumsum(rel_angle_nf);


%% plot stuff
%{
figure;
subplot(2, 1, 1);
time = (1:length(vd)) ./ fs;
plot(time, vd); hold on; plot(time, wiener_d); hold off;
title('Driver Velocity');
xlabel('time (s)');
ylabel('velocity (m/s)'); 
subplot(2, 1, 2);
time = (1:length(vp)) ./ fs;
plot(time, vp); hold on; plot(time, wiener_p); hold off;
title('Passenger Velocity');
xlabel('time (s)');
ylabel('velocity (m/s)'); 

% plot angle

figure;
time = (1:length(angle)) ./ fs;
plot(time, angle);
title('Vehicle Angle');
xlabel('time (s)');
ylabel('Angle (rad)');
%}

%% plot path

% Plot car path without noise no filtering 
v_ave = (vd + vp) / 2;
% v_ave = ones(length(angle), 1);
x = cos(angle) .* v_ave ./ fs;
y = sin(angle) .* v_ave ./ fs;
x_distance = cumsum(x);
y_distance = cumsum(y);

figure(1);
scatter(x_distance, y_distance, 'Marker','o');
hold on;


% Plot car path with noise no filtering 
v_ave = (vd_n + vp_n) / 2;
% v_ave = ones(length(angle), 1);
x = cos(angle_n) .* v_ave ./ fs;
y = sin(angle_n) .* v_ave ./ fs;
x_distance_n = cumsum(x);
y_distance_n = cumsum(y);

scatter(x_distance_n, y_distance_n, 'Marker','o');
title('Vehicle Path and Noise');
xlabel('Distance along x-axis (m)');
ylabel('Distance along y-axis (m)');
legend('Noise', 'No noise', 'Location', 'southeast'); 
% legend('No noise', 'Location', 'southeast'); 
% keyboard
% timedGraph(x_distance, y_distance);
 
% % Plot car path with noise and filtering
% v_ave = (wiener_d + wiener_p) / 2;
% % v_ave = ones(length(angle), 1);
% x = cos(angle_nf) .* v_ave ./ fs;
% y = sin(angle_nf) .* v_ave ./ fs;
% x_distance_nf = cumsum(x);
% y_distance_nf = cumsum(y);
% 
% 
% figure(3);
% scatter(x_distance_nf, y_distance_nf, 'Marker','o');
% title('Vehicle Path and Noise and filtering');
% xlabel('Distance along x-axis (m)');
% ylabel('Distance along y-axis (m)');
% hold on;

disp('SNR');
disp(snr);
disp('Average distance between points (m)'); % do rms instead
d = sum(sqrt((x_distance_n - x_distance).^2 + ((y_distance_n - y_distance).^2))) ./ length(x_distance);
disp(d);
disp('Error in Ending Location (m)');
e = sqrt((x_distance_n(end) - x_distance(end)).^2 + (y_distance_n(end) - y_distance(end)).^2);
dist = sum(sqrt((x_distance(2:end) - x_distance(1:end-1)).^2 + ((y_distance(2:end) - y_distance(1:end-1)).^2)));
disp(e);
disp('Distance Traveled (m)');
disp(dist);
disp('Percent Accuracy');
disp(abs(e - dist) / dist);


function timedGraph(x, y)
    figure; 
    timestep = 1;
    if (length(x) > 1000)
        timestep = int32(length(x) / 200);
    end
    scatter(x(1:1), y(1:1), 'Marker', 'o');
    title('Vehicle Path');
    xlabel('Distance along x-axis (m)');
    ylabel('Distance along y-axis (m)');
    xlim([-1, 30]);
    ylim([-1, 30]);
    pause(10);
    for i=1:timestep:length(x)
        scatter(x(1:i), y(1:i), 'Marker', 'o');
        xlim([-1, 30]);
        if (min(y) ~= max(y))
            ylim([-1, 30]);
        end
        title('Vehicle Path');
        xlabel('Distance along x-axis (m)');
        ylabel('Distance along y-axis (m)');
        pause(0.01);
    end
end
