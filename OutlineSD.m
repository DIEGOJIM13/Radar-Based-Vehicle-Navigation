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
theta = 30*pi/180; %z radar angle **********TEMP***************
d = 0.5; %dist between radar modules ************TEMP***********
lambda_d = physconst('LightSpeed')/fc_d;
lambda_p = physconst('LightSpeed')/fc_p;

%% Read in Radar Data
filename1 = strcat('4_9Data/test2/2.bin'); %open driver side txt file
filename2 = filename1; %open passanger side txt file

filename1 = strcat('tests/driver 4/0.bin'); %open driver side txt file
filename2 = strcat('tests/pass 4/0.bin'); %open driver side txt file

[a, b] = bin2txt(filename1, filename1);

% Driver
I1 = a(:, 2); % Get In phase data
Q1 = a(:, 3); % Get Quadrature data

% Passenger
I2 = b(:, 2); % Get In phase data
Q2 = b(:, 3); % Get Quadrature data



%% Recenter I and Q data and low pass

I1_t = I1 - mean(I1);
Q1_t = Q1 - mean(Q1);
I2_t = I2 - mean(I2);
Q2_t = Q2 - mean(Q2);

% fcut = 5000;
% [b,a] = butter(6,fcut/(fs/2));
% I1a = filter(b,a,I1_t);
% Q1a = filter(b,a,Q1_t);
% I2a = filter(b,a,I2_t);
% Q2a = filter(b,a,Q2_t);

%% convert to complex

dataCplx1 = complex(I1,Q1);
absArray1 = abs(dataCplx1);
phaseArray1 = angle(dataCplx1); %find phase - driver 

dataCplx2 = complex(I2,Q2);
absArray2 = abs(dataCplx2);
phaseArray2 = angle(dataCplx2);%find phase - passanger


%% convert to doppler frequency

fd = fs*(diff(phaseArray1)/(2*pi)); %driver doplar frequency

fp = fs*(diff(phaseArray2)/(2*pi)); %passanger doplar frequency

%% convert to velocity 
% radial velocity scaling by trig azimuth and elevation angle 

vd = (fd*c_air)/(2*fc_d); %find velocity driver
%     *cos(sqrt(2)/2)* sin(theta - 90) * -1; % scale for y and z and direction if car moving in x
 
vp = (fp*c_air)/(2*fc_p); %find velocity driver
%     *cos(sqrt(2)/2)* sin(theta - 90) * -1; % scale for y and z and direction if car moving in x

%% noise filtering

wiener_d = wiener2(vd, [1, 1000]);
wiener_p = wiener2(vp, [1, 1000]);

% wiener_d = vd;
% wiener_p = vp;

% Kalman Filtering - Currently not being used.
% s.A=1; s.Q=1e-7; s.H=1; s.R=1e-7; s.B=0; s.u=0; s.x=nan; s.P=nan;
% count = 0;
% for i = 1:1:length(velocity)
%     if (mod(i, int16(length(velocity)/100)) == 0)
%         count = count + 1;
%         disp(count);
%     end
%     s(end).z = velocity(i);
%     s(end + 1) = kalmanf(s(end));
% end
% 
% figure;
% subplot(2,1,1);
% plot(velocity);
% title('Original');
% 
% subplot(2,1,2);
% plot([s(2:end).x]);
% title('Kalman Filter');

%% calculate angle

omega = (wiener_d - wiener_p) ./ d;

rel_angle = omega * t;

angle = cumsum(rel_angle);


%% calculate distance 
% velocities between both and  

dd = wiener_d * t; % distance vd
dp = wiener_p * t;

dtd = cumsum(dd);
dtp = cumsum(dp);

td = sum(dd);
tp = sum(dp);
%disp(dt)

%% calculate distance 
% velocities between both and  

dd = vd*t;% distance vd
dp = vp*t;

dtd = cumsum(dd);
dtp = cumsum(dp);

td = sum(dd);
tp = sum(dp);
%disp(dt)

%% plot stuff

% plot normal velocity and wiener filter

% figure;
% subplot(2,2,1);
% plot(vd);
% title('Original');
% 
% subplot(2,2,2);
% [envHigh, envLow] = envelope(vd, 1000, 'peak'); envMean = (envHigh+envLow)/2;
% plot(vd); hold on; plot(envHigh);  hold on; plot(envLow);  hold on; plot(envMean);
% title('Original with envelope');
% 
% subplot(2,2,3);
% plot(wiener_d);
% title('Wiener filter');
% 
% subplot(2,2,4);
% [envHigh, envLow] = envelope(wiener_d, 1000, 'peak'); envMean = (envHigh+envLow)/2;
% plot(wiener_d); hold on; plot(envHigh);  hold on; plot(envLow);  hold on; plot(envMean);
% title('Wiener filter with envelope');

%plot I and Q data

% figure;
% time = (1:length(I1)) ./ fs;
% plot(time, I1);
% hold on;
% plot(time,Q1);
% title('I and Q');
% xlabel('time (s)');
% ylabel('I and Q value');

%plot velocity

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

%% plot path

% Plot car path
v_ave = (vd + vp) / 2;
% v_ave = ones(length(angle), 1);
x = cos(angle) .* v_ave ./ fs;
y = sin(angle) .* v_ave ./ fs;
x_distance = cumsum(x);
y_distance = cumsum(y);

figure;
scatter(x_distance, y_distance, 'Marker','o');
title('Vehicle Path');
xlabel('Distance along x-axis (m)');
ylabel('Distance along y-axis (m)');
hold on;

% Plot GPS points
% fileName = 'GPS Data/02-25-19-37-35.txt';
fileName = 'GPS Data/03-13-22-02-49.txt';
[Xs, Ys] = getGPSCoordinates(fileName, true);
scatter(Xs, Ys, 'Marker', 'x');
legend('Vehicle Location', 'GPS Data');
hold off;

% timedGraph(x_distance, y_distance);
% timedGraph(Xs, Ys);


%size = 5;
%xlims = [-size size]; ylims = [-size size];
%line([xlims nan 0 0],[0 0 nan ylims],'LineWidth',0.5, 'Color',[.2 .2 .2])
%axis square, grid on
%set(gca, 'XLim',xlims, 'YLim',ylims)
%title('Vehicle Path')
%xlabel('x distance (m)'), ylabel('y distance (m)')

%hold on

%plot(1,dtd, 'LineStyle','none', 'Marker','o', 'Color','m');


function [Xs, Ys] = getGPSCoordinates(fileName, usingStartandEnd)
    fileID = fopen(fileName,'r');
    lat = []; lon = [];
    started = false;
    while (true)
        line = fgetl(fileID);
        if (~ischar(line)) % Check that we haven't reached end of file.
            break;
        end
        if (length(line) <= 1) % Check that the line is not empty.
            continue;
        end
        if (contains(line, ' ')) % remove weird characters
           line = regexprep(line, ' ', '');
        end
        if (usingStartandEnd)
            if (strcmpi(line, 'Start')) % if we started data collection
                lat = []; lon = [];
                started = true;
            end
            if (strcmpi(line, 'End'))
                started = false;
            end
            if (~started) % if we have not started data collection.
               continue; 
            end
        end
        C = strsplit(line,',');
        if (length(C) < 3) % Check that the V or A character is there. 
            continue;
        end
        if (strcmpi(C{3}, 'V')) % Check that the data collection is valid
            continue;
        end
    %     disp(line);
        % convert degrees, minutes, and seconds to degrees
        d = strsplit(C{4},'.');
        len = d{1}(end-1:end);
        latD = str2double(d{1}(1:end-length(len)));
        latM = str2double(len);
        latS = str2double(d{2}) / 100;
        if (strcmp(C{5}, 'S'))
            latD = latD * -1;
        end
        d = strsplit(C{6},'.');
        len = d{1}(end-1:end);
        lonD = str2double(d{1}(1:end-length(len)));
        lonM = str2double(len);
        lonS = str2double(d{2}) / 100;
        if (strcmp(C{7}, 'S'))
            lonD = lonD * -1;
        end
        lat(end+1) = dms2degrees([latD, latM, latS]);
        lon(end+1) = dms2degrees([lonD, lonM, lonS]);
    end
    fclose(fileID);
    if (isempty(lat))
        return;
    end
    [Xs, Ys] = Spherical2AzimuthalEquidistant(lat, lon, 0, 0, 0, 0, 6371000 * pi);
    Xs = Xs - Xs(1);
    Ys = Ys - Ys(1);
end

function timedGraph(x, y)
    figure; 
    timestep = 1;
    if (length(x) > 1000)
        timestep = int32(length(x) / 100);
    end
    
    for i=1:timestep:length(x)
        scatter(x(1:i), y(1:i), 'Marker', '.');
        xlim([min(x), max(x)]);
        ylim([min(y), max(y)]);
        pause(1/100);
    end
end

