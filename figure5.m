clc
clear all
close all

% Simulated acoustic radiation force without glue for various Vpp (peak-to peak voltage)
Frad_1 = load("Frad_1vpp_corrected.txt");
Frad_2 = load("Frad_2vpp_corrected.txt");
Frad_3 = load("Frad_3vpp_corrected.txt");
Frad_4 = load("Frad_4vpp_corrected.txt");
Frad_5 = load("Frad_5vpp_corrected.txt");
Frad_6 = load("Frad_6vpp_corrected.txt");
Frad_7 = load("Frad_7vpp_corrected.txt");
Frad_8 = load("Frad_8vpp_corrected.txt");
Frad_9 = load("Frad_9vpp_corrected.txt");

% Simulated acoustic radiation force with glue for various Vpp
Frad_1_glue = load("Frad_1vpp_withglue.txt");
Frad_2_glue = load("Frad_2vpp_withglue.txt");
Frad_3_glue = load("Frad_3vpp_withglue.txt");
Frad_4_glue = load("Frad_4vpp_withglue.txt");
Frad_5_glue = load("Frad_5vpp_withglue.txt");
Frad_6_glue = load("Frad_6vpp_withglue.txt");
Frad_7_glue = load("Frad_7vpp_withglue.txt");
Frad_8_glue = load("Frad_8vpp_withglue.txt");
Frad_9_glue = load("Frad_9vpp_withglue.txt");

% Calculating trap stiffness for various Vpp without glue layer
[kt1,xdata1,ydata1,xfit1,yfit1] = comsol_stiffness(Frad_1,2.74);
[kt2,xdata2,ydata2,xfit2,yfit2] = comsol_stiffness(Frad_2,5.34);
[kt3,xdata3,ydata3,xfit3,yfit3] = comsol_stiffness(Frad_3,7.84);
[kt4,xdata4,ydata4,xfit4,yfit4] = comsol_stiffness(Frad_4,10.5);
[kt5,xdata5,ydata5,xfit5,yfit5] = comsol_stiffness(Frad_5,13.6);
[kt6,xdata6,ydata6,xfit6,yfit6] = comsol_stiffness(Frad_6,16.2);
[kt7,xdata7,ydata7,xfit7,yfit7] = comsol_stiffness(Frad_7,18.7);
[kt8,xdata8,ydata8,xfit8,yfit8] = comsol_stiffness(Frad_8,21.2);
[kt9,xdata9,ydata9,xfit9,yfit9] = comsol_stiffness(Frad_9,23.7);

% Calculating trap stiffness for various Vpp with glue layer
[kt1_glue,xdata1_glue,ydata1_glue,xfit1_glue,yfit1_glue]  = comsol_stiffness(Frad_1_glue,2.74);
[kt2_glue,xdata2_glue,ydata2_glue,xfit2_glue,yfit2_glue] = comsol_stiffness(Frad_2_glue,2.74);
[kt3_glue,xdata3_glue,ydata3_glue,xfit3_glue,yfit3_glue] = comsol_stiffness(Frad_3_glue,2.74);
[kt4_glue,xdata4_glue,ydata4_glue,xfit4_glue,yfit4_glue] = comsol_stiffness(Frad_4_glue,2.74);
[kt5_glue,xdata5_glue,ydata5_glue,xfit5_glue,yfit5_glue] = comsol_stiffness(Frad_5_glue,2.74);
[kt6_glue,xdata6_glue,ydata6_glue,xfit6_glue,yfit6_glue] = comsol_stiffness(Frad_6_glue,2.74);
[kt7_glue,xdata7_glue,ydata7_glue,xfit7_glue,yfit7_glue] = comsol_stiffness(Frad_7_glue,2.74);
[kt8_glue,xdata8_glue,ydata8_glue,xfit8_glue,yfit8_glue] = comsol_stiffness(Frad_8_glue,2.74);
[kt9_glue,xdata9_glue,ydata9_glue,xfit9_glue,yfit9_glue] = comsol_stiffness(Frad_9_glue,2.74);

% storing the trap stiffnesses without glue layer
kt = [kt1 kt2 kt3 kt4 kt5 kt6 kt7 kt8 kt9];

% storing the trap stiffnesses with glue layer
kt_glue = [kt1_glue kt2_glue kt3_glue kt4_glue kt5_glue kt6_glue kt7_glue kt8_glue kt9_glue]; 

% store the applied Vpp values
vpp = [2.74 5.34 7.84 10.5 13.6 16.2 18.7 21.2 23.7];


% plot the sinusoidal fit to simulated radiation force
plt = Plot(xfit9, yfit9, xfit9_glue,yfit9_glue);
hold on
plot(xdata9,ydata9,'o','LineWidth',1.5)
plot(xdata9_glue,ydata9_glue,'o','LineWidth',1.5)
plt.BoxDim = [4, 3]; %[width, height]
plt.XLabel = 'z (\mum)'; % xlabel
plt.YLabel = 'Frad (pN)'; %ylabel
plt.XMinorTick = 'off';
plt.YMinorTick = 'off'; 
hold off


% plot the stiffness vs. vpp curve
plt = Plot(vpp,kt,vpp,kt_glue);
hold on
plot(vpp,kt,'o','LineWidth',1.5)
plot(vpp,kt_glue,'o','LineWidth',1.5)
plt.BoxDim = [4, 3]; %[width, height]
plt.XLabel = 'V_{pp} (V)'
plt.YLabel = 'Trap stiffness (pN/\mum)'; %ylabel
plt.XMinorTick = 'off'; % 'on' or 'off'
plt.YMinorTick = 'off'; % 'on' or 'off'
hold off

% Function that calculates the trapping stiffness using simulated data
function [kt,xdata,ydata,xfit,yfit] = comsol_stiffness(data,n)

% storing x and y data for plot
x = data(:,2); % position in the channel
y = data(:,3); % Acoustic radiation force 

% setting up initial values to fit sinusoidal to acoustic radiation force
yu = max(y);
yl = min(y);
yr = (yu-yl);                               % Range of ‘y’
yz = y-yu+(yr/2);
zx = x(yz .* circshift(yz,[0 1]) <= 0);     % Find zero-crossings
per = 100;                     % Estimate period
ym = mean(y);                               % Estimate offset

% fitting sinusoidal to acoustic radiation force
fit = @(b,x)  b(1).*(sin(2*pi*x./(1/per) + 2*pi/b(2))) + b(3);    % Function to fit
fcn = @(b) sum((fit(b,x) - y).^2);                              % Least-Squares cost function
s = fminsearch(fcn, [yr/2;-1;  ym])                       % Minimise Least-Squares
xp = linspace(min(x),max(x));
yfit = fit(s,xp)*1e12; % convert fitted data from N to pN
xfit = xp+100; % setting bottom of the channel to be x=0 micrometer 

% converting simulation data to correct reference and units
ydata = data(:,3)*1e12; % convert simulation data from N to pN
xdata = data(:,2)+100; % setting bottom of the channel to be x=0 micrometer 

% find zero crossing
y_val = fit(s,xp);
y_val_cross = y_val(1:end-1).*y_val(2:end);
[val,ind] = min(y_val_cross);

% calculate trap stiffness
stiffness = abs((y_val(ind+1) - y_val(ind-1))/((xp(ind+1)-xp(ind-1))*1e-6)); % N/m
kt = stiffness*1e6; % convert trap stiffness to pN/micrometer
end

