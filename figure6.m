% MATLAB code that calculate and plot the calibration curve

clc
clear all
close all

% constant image processing parameters
level = 0.30;
width = 50;
height = 50;


% track beads in in focus calibration image
file = 'Calibration sample\1ms dial led 30 2\calibration_30dialled_64analoggain_exposure1ms_470nmwavelength_dialforward 0.tif';

% run calibration
[resnorms, output] = calib(file,width,height,level);

% setting the objective moving up in axial direction as positive
z_pos = -output(:,2);
fwhm = output(:,1);

% removing unnecessary data
[ind1,~] = find(z_pos>5);
z_pos(ind1) = [];
[ind2,~] = find(z_pos<-5);
z_pos(ind2) = [];
fwhm(ind1) = [];
fwhm(ind2) = [];


% calculate the median radii for each axial position
min_z = min(z_pos);
max_z = max(z_pos);
z_med =  z_pos - min_z + 1;
median_fwhm = accumarray(z_med,fwhm,[],@median);
std_fwhm = accumarray(z_med,fwhm,[],@std);


% removing zeros from median radii, since not all position got a value and
% get assigned zero 
median_fwhm = median_fwhm(median_fwhm~=0)';
std_fwhm = std_fwhm(std_fwhm~=0)';

% create a vector with unique axial positions
z_val = sort(unique(z_pos))';

% fit a fourth order polynomial to the calibration data 
z_fit = -6:0.01:6;
fit_eq = ' p1*x^4 + p2*x^2 + p3';
p1_es = 1;
p2_es = 1;
p3_es = min(median_fwhm);
startPoints = [p1_es p2_es p3_es];
tf = excludedata(z_val,median_fwhm,'domain',[-5 5]);
f = fit(z_val',median_fwhm',fit_eq,'Start', startPoints,'Exclude', tf,'Lower',[0,0,0]);

% Evaluate the fit function for finer axial position
p1 = f.p1;
p2 = f.p2;
p3 = f.p3;
rad_fit = p1*z_fit.^4 + p2*z_fit.^2 + p3;

% calculate R^2 of the fit
y_fit = p1*z_val.^4 + p2*z_val.^2 + p3;
y_val = median_fwhm;
SStot = sum((y_val-mean(y_val)).^2);                    % Total Sum-Of-Squares
SSres = sum((y_val-y_fit).^2);                       % Residual Sum-Of-Squares
Rsq = 1-SSres/SStot;                            % R^2

% plotting calibration curve
plt = Plot(z_fit,rad_fit);
hold on
boxchart(z_pos,fwhm)
plt.BoxDim = [4, 3];
%plt.Title = 'Calibration curve'; % plot title
% plt.XLabel = 'z (\mum) distance away from focal plane'; % xlabel
plt.XLabel = 'z (\mum)'; % xlabel
plt.YLabel = '\sigma_x (\mum)'; %ylabel
%legend('Fit','data')
%grid on


% functions
function bead = bead_tracking(filename,level)
A = imread(filename);

% remove scale bar
A(1418:1445,1266:1310,:) = 0; % calibration sample

% turn the calibration image to gray image
A = rgb2gray(A);

% show gray image
figure
imshow(A)
title ('gray image of the beads')

% convert to black and white image
A = im2bw(A,level);

% show image of black and white image
% figure
% imshow(A)
% title('Black and white image')

%imtool(A)
% 16 pixels equals 2 micrometers from imtool picture so 1.5 micrometer is
% equal to:
ratio = 2/16;

% extra filter to filter out right size particles
figure
imshow(A)
title('Detected beads')

% we let the size vary by 10 percent so the range is
range = [1 10];

[centers, radii, metric] = imfindcircles(A,range);

%Draw circels around detected beads
viscircles(centers, radii,'EdgeColor','b');

% save center and radii of the detected beads
bead = [centers, radii];
end

function [width_error, resnorm, fwhm,fit_param] = fit2DGaussian(J,bead,width,height,xdata, lb, ub)

% Initial fit parameters
amp = max(J,[],'all');
wx = bead(1,3);
wy = wx;
x_center = bead(1,1);
y_center = bead(1,2);
x0 = [amp, x_center, wx , y_center, wy];

% fitting
[x,resnorm,residual] = lsqcurvefit(@D2GaussFunction,x0,xdata,J,lb,ub);
fit_param = x;
sig_x = x(3);
sig_y = x(5);
width_error = abs(sig_x-sig_y);
fwhm = sig_x;
end



function [resnorms,output] = calib(filename,width,height,level)

bead = bead_tracking(filename,level);
[beadRows,beadCols] = size(bead);


ratio = 2/16; % convert pixels to micrometers

% setting up matrices to store data
data = [];
output = [];
resnorm_mat = []; 
resnorms = [];
temp = 0;

% create meshgrid to plot 2D Gaussian
[X,Y] = meshgrid(-width/2:height/2);
xdata = zeros(size(X,1),size(Y,2),2);
xdata(:,:,1) = X;
xdata(:,:,2) = Y;

% setting upper and lower for bound on Guassian fit
lb = [0,-width/2,0,-height/2,0];
ub = [realmax('double'),width/2,(width/2)^2,height/2,(width/2)^2];

% loop through all the beads
for j = 1:beadRows
    temp = temp + 1;

    % calculate the cropped image size
    xmin = round(bead(j,1)-0.5*width);
    xmax = round(bead(j,1)+0.5*width);
    ymin = round(bead(j,2)-0.5*height);
    ymax = round(bead(j,2)+0.5*height);

    % loop through the z stack
    for n = 11:-1:-7

        str = 'Calibration sample\1ms dial led 30 2\calibration_30dialled_64analoggain_exposure1ms_470nmwavelength_dialforward ';
        filename = append(str,num2str(n),'.tif' );

        % converting RGB to gray image
        grayImage = rgb2gray(imread(filename));

        % subtracting the mean
        grayImage = grayImage - mean(grayImage(:));

        % calculate the size of the frame
        [rows, columns, numberOfColorChannels] = size(grayImage);

        % check if cropping is within the bounds of the original image
        if ymax<=size(grayImage,1) && ymax>=1 && ymin>=1 && ymin<=size(grayImage,1)&& xmax<=size(grayImage,2) && xmax>=1 && xmin>=1 && xmin<=size(grayImage,2)

            % crop image
            J_gray = grayImage(ymin:ymax,xmin:xmax);

            % turn image format to double
            J = double(J_gray);

            % normalize intensity profile
            J = J/max(J,[],'all');

            if sum(J(:))~=0

                % Fit bead with 2D Gaussian function
                [width_error, resnorm, fwhm,fit_param] = fit2DGaussian(J,bead(j,:),width,height,xdata, lb, ub);

                % filter out bad data
                if resnorm<32 && width_error < 1 %&& abs(max(D2GaussFunction(fit_param,xdata),[],"all")-max(J(:)))<0.35

                    % save fwhm and z position of bead
                    fwhm = fit_param(3)*ratio;
                    %z_pos = 2*n;
                    z_pos = n;
                    bead_data = [fwhm,z_pos];
                    data = [data;bead_data];
                    resnorm_mat = [resnorm_mat;resnorm];
                    
%                     J_mat = zeros(size(J,1),size(J,2),19);
%                     if j ==1
%                     gsFilename = sprintf('%d distance away focal plane.jpg', -1*(n+2));
%                     imwrite(J_gray,gsFilename);
% %                     figure
% %                     imshow(J_gray)
%                     end

%                     % plotting
%                     figure
%                     C = del2(J);
%                     mesh(X,Y,J,C) %plot data
%                     mesh(X,Y,J,'EdgeColor', 'r')
%                     %mesh(X,Y,J) %plot data
%                     hold on
%                     surface(X,Y,D2GaussFunction(fit_param,xdata),'EdgeColor','b') %plot fit
%                     axis([-width/2-0.5 width/2+0.5 -height/2-0.5 height/2+0.5])
%                     legend('Data','fit')
%                     str = 'Resnorm of %d';
%                     str = sprintf(str,resnorm);
%                     title(str)
% 
%                     alpha(0.2)
%                     hold off

                end
            end
        end
    end

    % Take the smallest bead size as on the focal plane

    if ~isempty(data) && length(data(:,1))>3
        z = data(:,2);
        fwhm_val = data(:,1);

        %[min_fwhm,ind] = min(data_filt(:,1));
        [min_fwhm,ind] = min(fwhm_val);
        z_cor = z-z(ind);

        % save bead data
 
        out = [fwhm_val,z_cor];
        output = [output;out];
        resnorms = [resnorms;resnorm_mat];
    
    end

    % empty data
    data = [];
    resnorm_mat = [];
end
end