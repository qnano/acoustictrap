% MATLAB code that derive the variance in z (axial coordinate) of the
% microbead and calculates the corresponding trap stiffness.

% Note: This code uses the saved output
% data calculated in bead_experiments_fwhm.m   

clc
clear all
close all

% loading saved sigma_x derived with bead_experiments_fwhm.m   

output_3vpp = load('experimental data/dschip_3vpp.mat');
[kt_3vpp_max kt_3vpp  kt_3vpp_min error_3vpp] = stiffness_calc(output_3vpp.output,7.84)


output_4vpp = load('experimental data/dschip_4vpp.mat');
[kt_4vpp_max kt_4vpp kt_4vpp_min error_4vpp] = stiffness_calc(output_4vpp.output,10.5)


output_5vpp = load('experimental data/dschip_5vpp.mat');
[kt_5vpp_max kt_5vpp kt_5vpp_min error_5vpp] = stiffness_calc(output_5vpp.output,13.6)


output_6vpp = load('experimental data/dschip_6vpp.mat');
[kt_6vpp_max kt_6vpp kt_6vpp_min error_6vpp] = stiffness_calc(output_6vpp.output,16.2)

output_7vpp = load('experimental data/dschip_7vpp.mat');
[kt_7vpp_max kt_7vpp kt_7vpp_min error_7vpp] = stiffness_calc(output_7vpp.output,18.7)

output_8vpp = load('experimental data/dschip_8vpp.mat');
[kt_8vpp_max kt_8vpp kt_8vpp_min error_8vpp] = stiffness_calc(output_8vpp.output,21.2)


output_9vpp = load('experimental data/dschip_9vpp.mat');
[kt_9vpp_max kt_9vpp kt_9vpp_min error_9vpp] = stiffness_calc(output_9vpp.output,23.7);

% setting up stiffness vectors and peak-to-peak voltage vector
kt =  [kt_3vpp kt_4vpp kt_5vpp kt_6vpp kt_7vpp kt_8vpp kt_9vpp];
kt_max =  [kt_3vpp_max kt_4vpp_max kt_5vpp_max kt_6vpp_max kt_7vpp_max kt_8vpp_max kt_9vpp_max];
kt_min =  [kt_3vpp_min kt_4vpp_min  kt_5vpp_min  kt_6vpp_min  kt_7vpp_min  kt_8vpp_min  kt_9vpp_min];
vpp = [7.84 10.5 13.6 16.2 18.7 21.2 23.7];

% loading trap stiffness data derived with COMSOL
kt_comsol = load('kt_comsol_corrected.mat'); % trap stiffness without glue
kt_comsol_glue = load('kt_comsol_glue.mat'); % trap stiffness with glue

% create publication quality plots
xdata = vpp;
ydata = kt;

% fitting the experimental obtained stiffness
p = polyfit(xdata,ydata,2);
xfit = linspace(min(xdata),max(xdata),100);
yfit = polyval(p,xfit);

% setting up vector to contain the uncertainty in trap stiffness and
% voltage
error_kt = [error_3vpp error_4vpp error_5vpp error_6vpp error_7vpp error_8vpp error_9vpp];
error_Vpp = vpp*0.01; % 1 percent of uncertainty in Vpp values 

% setting up vector to contain the simulated trap stiffness with and
% without glue layer
ycomsol = kt_comsol.kt(3:end);
ycomsol_glue = kt_comsol_glue.kt_glue(3:end)

% plot trap stiffnesses from model and experiments
plt = Plot(xfit, yfit);
hold on
plot(xdata,ycomsol,'--','LineWidth',1.5);
plot(xdata,ycomsol_glue,'b--','LineWidth',1.5);
e = errorbar(xdata, ydata,error_kt,error_kt,error_Vpp,error_Vpp,'LineWidth',1,'Color' ,'k');
e.LineStyle = 'none';


% plot settings
plt.BoxDim = [4, 3]; %[width, height]
%plt.Title = 'Experimental vs. Simulated stiffness'; % plot title
plt.XLabel = 'V_{pp} (V)'; % xlabel
plt.YLabel = 'Trap stiffness (pN/\mum)'; %ylabel
%ax.GridAlpha = 0.1;
%grid on
%legend('Experiment','Simulation','Simulation with glue')
hold off

% code that calculates the trap stiffness from sigma_x data
function [kt_max kt kt_min error] = stiffness_calc(output,n)

% setting output equal
fwhm_beads = output;

%choose number of bins based on Sturge's Rule
nbins = round(1+3.322*log(length(output)));

% getting the fwhm value with highest frequency in histogram
[N,edges] = histcounts(fwhm_beads,nbins);
[maxN,ind] = max(N);
M = edges(ind+1)
fwhm_beads(fwhm_beads<M) = [];


% loading calibration curve parameters
load('p_val_augustlast.mat')

% adjust the origin of the calibration curve to reflect the correct focal
% plane
p3 = M; 

z_val = [];

for i = 1:length(fwhm_beads)

        % find the z value that corresponds to the fwhm value
        val = fwhm_beads(i);
        pol = [ p1  0 p2  0 p3-val];
        r = roots(pol);
        r = r(imag(r)==0);
        r = r(r>=0); 
        z= r;
    z_val = [z_val;z];
end


% convert output data from pixels to micrometers and get rid of too many
% decimals
rounding = 100; 

% determine z value with highest frequency and set that to focal plane
% needs to be done due to too many decimals in MATLAB that makes noise data 
z_val1 = round(z_val*rounding)/rounding;
M = mode(z_val1);
z_val(z_val<M) = [];

% set z value with highest frequency equal to zero
z_val = z_val - M;

% setting up symmetric z values
z_val = [-z_val(z_val~=0)-max(-z_val(z_val~=0));z_val];

% calculate the required number of bins again using Sturge's Rule for this
nbins = round(1+3.322*log(length(z_val)));

% calculate the frequency and bin edges of histogram
[N,edges] = histcounts(z_val,nbins);

% calculate the midpoint of each histogram bin 
x = (edges(1:end-1)+edges(2:end))/2;

% remove midpoint where there is zero frequency
x = x(N~=0);

% remove all the zero frequencies of the histogram
N = N(N~=0);

% center around the maximum frequency  
[~,ind]  = max(N);
edges = edges-x(ind); % adjust center of the bin with highest frequency to be zero

[N,edges] = histcounts(z_val,edges) % recalculate counts again

% remove all the zero frequencies of the histogram
N = N(N~=0);

% calculate the midpoint of each histogram bin again 
x = (edges(1:end-1)+edges(2:end))/2;

% remove midpoint where there is zero frequency
x = x(N~=0);

% center around the maximum frequency  
[~,ind]  = max(N);
edges = edges-x(ind); % adjust center of the bin with highest frequency to be zero
x = x-x(ind);
% [N,~] = histcounts(z_val,edges)

% fitting Gaussian to histogram
norm_fun = @(sd,x) max(N)*exp(-(x).^2/(2*sd^2));
[sol, f] = lsqcurvefit(@(p, x) norm_fun(p(1), x), [rand(1,1)], x,N );


% get the standard deviation in z position
sd = sol;

if n == 18.7
%  plot z values as histogram with fit
x_p = min(x):0.01:max(x);
plt = Plot(x_p, norm_fun(sd,x_p));
hold on
histogram(z_val,edges)
%plot(x,N,'o')
plt.BoxDim = [4, 3]; 
plt.XLabel = 'z (\mum)'; % xlabel
plt.YLabel = 'Frequency'; %ylabel
str = 'Gaussian fit histogram of z values for %d Vpp';
str = sprintf(str,n);
plt.Title = str;
%ax.GridAlpha = 0.1;
%grid on
legend('Fit','Experiment')
end

% Boltzman constant 
kB = 1.380649e-23 ;  %1.380649 × 10-23 m2 kg s-2 K-1

% Assumed temperature in channel
T = 293.15; % 20 degrees celsius in Kelvin

% calculate trap stiffness with equipartition theory 
trapp_stiffness = 2*kB*T/(sd*1*10^(-6))^2; % trapp stiffness in N/m
kt = trapp_stiffness*1e6; % trapp stiffness in pN/micrometer

% calculate the uncertainty in experimental trap stiffness
delta_T = 2; % error in temperature
delta_wx = 0.06; % error in wx values (FWHM) (should be the mean error in z range of -5 to 5)
wx = 0.39; % wx value at z = zero (should be the mean in z range of -5 to 5)
delta_b = 0.022; % error in bead size
b = 1.377;

frac_error_var = delta_wx/wx + delta_b/b;
frac_error_T = delta_T/T;
error = (frac_error_var + frac_error_T)*kt;
kt_max = kt + error; 
kt_min = kt - error;
end 

















