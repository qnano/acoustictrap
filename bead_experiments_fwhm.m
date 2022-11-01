% This code calculates the sigma_x value in bead trapping videos using 2D
% Gaussian fit to the intensity profile of each bead in the user defined
% ROI

% Note: You need to provide the correct path to the trapping video and
% tracking software to make this code work correctly 

clc
clear all
close all
imtool close all

% track beads in bead trapping video (add path to trap videos folder)
file = '1-8 measurements\dschiplast1-8-2022\last_mearumentdschip9vpp.mp4' 

%file = 'E:\Thesis MATLAB\Chip stiffness_best\MATLAB codes\1-8 measurements\dschiplast1-8-2022\last_mearumentdschip9vpp.mp4';
% image cropping constants
width = 20;
height = 20;
numFrames = 0;

% image processing values
peak_threshold = 20; % find peaks above certain value
average_diam = 5;
threshold = 1;
diam_est = 5; %(and 16) estimate of diameter of beads
freq_cutoff = 1; % spatial wavelength cutt off in pixels almost always 1
window_size = 7;
nFrame = 1; 

% value to choose between analysis on a single frame for val == 1 or all
% video frames vall ~= 1
val = 0;

% read trapping video
v = VideoReader(file);

% calculating number of frames in 30 seconds video
while v.CurrentTime <= 30
    readFrame(v);
    numFrames = numFrames + 1;
end



% add path to image processing software
addpath 'E:\Thesis MATLAB\Chip stiffness_best\MATLAB codes\tracking software'; % laptop
%addpath 'D:\Thesis MATLAB\Chip stiffness_best\tracking software\'; % desktop

tic
output = trap_stiffness_fun(file,width,height,peak_threshold,average_diam,diam_est,freq_cutoff,val,nFrame,window_size,threshold);
toc

% plot histogram of sigma_x values
figure
histogram(output)



% 2D Gaussian function to calculate to calculate sigma_x/fwhm from the trap videos
function [width_error, resnorm, fwhm,fit_param] = fit2DGaussian(J,bead,xdata, lb, ub)
amp = max(J,[],'all');
wx = bead(1,3);
wy = wx;
x_center = bead(1,1);
y_center = bead(1,2);
bg = median(J); % background


% Initial fit parameters
x0 = [amp, x_center, wx , y_center, wy, bg];

% fitting
[x,resnorm,residual] = lsqcurvefit(@D2GaussFunction,x0,xdata,J,lb,ub);
fit_param = x;
sig_x = x(3);
sig_y = x(5);
width_error = abs(sig_x-sig_y);
fwhm = sig_x;
end

% function to perform 2D Gaussian fit on each bead in the ROI
function output = trap_stiffness_fun(filename,width,height,peak_threshold,average_diam,diam_est,freq_cutoff,val,nFrame,window_size,threshold)

% setting up matrices to store data
data = [];
output = [];
temp = 0;
numFrames = 0;
ratio = 2/5;

% create meshgrid to plot 2D Gaussian
[X,Y] = meshgrid(-width/2:height/2);
xdata = zeros(size(X,1),size(Y,2),2);
xdata(:,:,1) = X;
xdata(:,:,2) = Y;

% setting upper and lower for bound on Guassian fit
lb = [0,-width/2,0,-height/2,0];
ub = [realmax('double'),width/2,(width/2)^2,height/2,(width/2)^2];

% read trapping video
v = VideoReader(filename);

% calculating number of frames in 30 seconds video
while v.CurrentTime <= 30
    readFrame(v);
    numFrames = numFrames + 1;
end

% Get first frame for user to crop
frame = read(v,nFrame);


% analyse beads with imtool
%imtool(frame)

% let user crop the frame
figure
imshow(frame)
message = sprintf('Left click and hold to begin drawing.\nClose the image to see the cropped image');
uiwait(msgbox(message));

% hFH = imfreehand();
h_rect = imrect();

pos_rect = h_rect.getPosition();

% Round off so the coordinates can be used as indices
pos_rect = round(pos_rect);

% analyse bead tracking and Gaussian fitting
if val == 1

    frame = read(v,nFrame);

    % converting RGB to gray image
    grayImage = rgb2gray(frame);

    % crop grayImage and frame as specified by user for all other frames
    grayImage = grayImage(pos_rect(2) + (0:pos_rect(4)), pos_rect(1) + (0:pos_rect(3)));
    frame = frame(pos_rect(2) + (0:pos_rect(4)), pos_rect(1) + (0:pos_rect(3)));


    a = double(grayImage);
    b = bpass(a,freq_cutoff,diam_est,threshold);
    pk = pkfnd(b,peak_threshold,average_diam);
    bead = cntrd(b,pk,window_size);

    figure
    imshow(frame)

    figure
    imshow(frame);
    title('Cropped frame')
    hold on
    viscircles(bead(:,1:2), bead(:,3)/2000,'EdgeColor','b');
    hold off


    [beadRows,beadCols] = size(bead);


    % subtracting the mean
    grayImage = grayImage - mean(grayImage(:));

    % counter to store data
    temp = temp + 1;

    % loop through all the beads
    for j = 1:beadRows


        % calculate the cropped image size
        xmin = round(bead(j,1)-0.5*width);
        xmax = round(bead(j,1)+0.5*width);
        ymin = round(bead(j,2)-0.5*height);
        ymax = round(bead(j,2)+0.5*height);


        % calculate the size of the frame
       [rows, columns, numberOfColorChannels] = size(grayImage);

        % check if cropping is within the bounds of the original image
        if ymax<=size(grayImage,1) && ymax>=1 && ymin>=1 && ymin<=size(grayImage,1)&& xmax<=size(grayImage,2) && xmax>=1 && xmin>=1 && xmin<=size(grayImage,2)


            % crop image
            J = grayImage(ymin:ymax,xmin:xmax);

            figure
            imshow(J)
            title('Cropped image around bead')

            % convert cropped image to double
            J = double(J);

            % normalize cropped image
            J = J/max(J,[],'all');



            if sum(J(:))~=0

                % Fit bead with 2D Gaussian function
                if isequal(size(J,1), size(xdata,1)) && isequal(size(J,2), size(xdata,2))
                    [width_error, resnorm, fwhm,fit_param] = fit2DGaussian(J,bead(j,:),xdata, lb, ub);


                    % filter out bad data
                    if resnorm<32 && width_error < 1 %&& abs(max(D2GaussFunction(fit_param,xdata),[],"all")-max(J(:)))<0.35

                        % save fwhm and z position of bead
                        fwhm = fit_param(3)*ratio;
                        data = [data;fwhm];


                        % plotting
                        figure
                        C = del2(J);
                        mesh(X,Y,J,C) %plot data
                        mesh(X,Y,J,'EdgeColor', 'r')
                        %mesh(X,Y,J) %plot data
                        hold on
                        surface(X,Y,D2GaussFunction(fit_param,xdata),'EdgeColor','b') %plot fit
                        axis([-width/2-0.5 width/2+0.5 -height/2-0.5 height/2+0.5])
                        legend('Data','fit')
                        alpha(0.2)
                        hold off

                    end
                else
                    continue
                end
            end
        end
    end



    % store data
    output = [output;data];

    % empty data
    data = [];

else
    % loop through video frames
    for i = 1:numFrames

        frame = read(v,i);


        % converting RGB to gray image
        grayImage = rgb2gray(frame);

        % crop grayImage and frame as specified by user for all other frames
        grayImage = grayImage(pos_rect(2) + (0:pos_rect(4)), pos_rect(1) + (0:pos_rect(3)));
        %  frame = frame(pos_rect(2) + (0:pos_rect(4)), pos_rect(1) + (0:pos_rect(3)));


        a = double(grayImage);
        b = bpass(a,freq_cutoff,diam_est,threshold);
        pk = pkfnd(b,peak_threshold,average_diam);
        bead = cntrd(b,pk,window_size);

        %     figure
        %     imshow(grayImage);
        %     title('Cropped frame')
        %     hold on
        %     viscircles(bead(:,1:2), bead(:,3)/1000,'EdgeColor','b');
        %     hold off



        [beadRows,beadCols] = size(bead);


        % subtracting the mean
        grayImage = grayImage - mean(grayImage(:));

        % counter to store data
        temp = temp + 1;

        % loop through all the beads
        parfor j = 1:beadRows


            % calculate the cropped image size
            xmin = round(bead(j,1)-0.5*width);
            xmax = round(bead(j,1)+0.5*width);
            ymin = round(bead(j,2)-0.5*height);
            ymax = round(bead(j,2)+0.5*height);


            % calculate the size of the frame
            [rows, columns, numberOfColorChannels] = size(grayImage);

            % check if cropping is within the bounds of the original image
            if ymax<=size(grayImage,1) && ymax>=1 && ymin>=1 && ymin<=size(grayImage,1)&& xmax<=size(grayImage,2) && xmax>=1 && xmin>=1 && xmin<=size(grayImage,2)


                % crop image
                J = grayImage(ymin:ymax,xmin:xmax);


                % convert cropped image to double
                J = double(J);

                % normalize cropped image
                J = J/max(J,[],'all');


                if sum(J(:))~=0

                    % Fit bead with 2D Gaussian function
                    if isequal(size(J,1), size(xdata,1)) && isequal(size(J,2), size(xdata,2))
                        [width_error, resnorm, fwhm,fit_param] = fit2DGaussian(J,bead(j,:),xdata, lb, ub);


                        % filter out bad data
                        if resnorm<32 && width_error < 1 %&& abs(max(D2GaussFunction(fit_param,xdata),[],"all")-max(J(:)))<0.35

                            % save fwhm and z position of bead
                            fwhm = fit_param(3)*ratio;
                            data = [data;fwhm];


                            % plotting
                            %                         figure
                            %                         C = del2(J);
                            %                         mesh(X,Y,J,C) %plot data
                            %                         mesh(X,Y,J,'EdgeColor', 'r')
                            %                         %mesh(X,Y,J) %plot data
                            %                         hold on
                            %                         surface(X,Y,D2GaussFunction(fit_param,xdata),'EdgeColor','b') %plot fit
                            %                         axis([-width/2-0.5 width/2+0.5 -height/2-0.5 height/2+0.5])
                            %                         legend('Data','fit')
                            %                         alpha(0.2)
                            %                         hold off

                        end
                    else
                        continue
                    end
                end
            end
        end


        % store data
        output = [output;data];

        % empty data
        data = [];

    end
end

end








