%% Dov Pupil Code
%This script will find the pupil diameter and position in a video of the
%mouse eye. Modified from a version written by Matt McGinley.

%required functions:
%fit_ellipse
%ellipse
%inpaint
%Image Acquisition Toolbox and all that good stuff

%Basic outline to script
%1) Load video of eye (3D matrix) and crop (before for loop).
%2) Select IR light reflection and fill in by extrapolation (before for loop). REMOVED
%3) Find dark/pupil pixels with user value threshold. Estimate pupil center.
%and radius. Delete pixels which are too far from center to be pupil.
%4) Canny edge detection on original image.
%5) Delete edge pixels which are too far from center and dark pixels. Also
%delete pixels which are not continguous with many other edge pixels.
%6) Fit ellipse to edge pixels
%7) Fill in any values when pupil could not be fit (e.g. eye closing).

%% Load video and crop

%video is a r x c x t  uint8 matrix where r is height, c is width, and t is video frame
%videoTS timestamps of each video frame
load('video.mat')
load('videoTS.mat')

%get rectangle for cropping
%drag rectangle over just eye region
figure
uiwait(msgbox('Drag box over eye -> Right click -> Crop Image'));
[~,rectEye] = imcrop(video(:,:,1));
close all
%crop every image in movie
video = video(round(rectEye(2)):round(rectEye(2)+rectEye(4)),round(rectEye(1)):round(rectEye(1)+rectEye(3)),:);
disp('video cropped');

% %% Remove IR light reflection as best we can (optional)
% %Test random frames before extrapolating on all frames and saving
%
% %manually select outline of reflection
% im = video(:,:,1);
% maskLight = roipoly(im);
% close all
%
% for iTest = 1:10
%     randFrame = ceil(rand*size(video,3));
%     disp(['random frame #' num2str(randFrame)])
%     randImg = video(:,:,randFrame);
%     J = regionfill(randImg,maskLight);
%     K = edge(randImg,'canny');
%     L = edge(J,'canny');
%     subplot(2,2,1)
%     imshow(randImg)
%     subplot(2,2,3)
%     imshow(J)
%     subplot(2,2,2)
%     imshow(K)
%     subplot(2,2,4)
%     imshow(L)
%     pause
% end

%% Select eye region of interest to improve quality

%average image of pupil
picEyeAvg = mean(video,3);
figure
uiwait(msgbox('Select region of interest -> Right click -> Create Mask '));
maskEye = roipoly(picEyeAvg./255);
close all

%% Choose threshold for pupil value

threshDark = 30;

figure
uiwait(msgbox('Choose threshDark value so that most pixels of pupil are showing. Press ESC to progress.'));
for iTest = 1:15
    randFrame = ceil(rand*size(video,3));
    disp(['random frame #' num2str(randFrame)])
    randImg = video(:,:,randFrame);
    randImg_threshed = randImg<threshDark;
    imshow(randImg_threshed.*maskEye)
    pause
end
close all

%% Begin

%eye may be closed if very few dark pixels found
threshNdarkPix = 200; %default 200
%fudge factor for estimated radius
radiusFudge = 1.4;  %default 1.4
%edge detection threshold fudge factor
edgeFudge = 1;   %default 1 (lower means more edges will be found)
%edge pixels must be at least this close to dark pixels
threshPupilDist = 3;  %default 3
%edge pixels must be continguous with at least this many other pixels
minEdgeSize = 15;     %default 15
%plot and pause after each frame?
plot_each_frame=false;

pupilD = nan(1,size(video,3));
pupilXY = nan(2,size(video,3));
%videoAnnotated = zeros(388,538,3,5000,'uint8');

tic
for iFrame = 1:size(video,3)
    %for ploting pupil edges
    if plot_each_frame==true
        row3 = [];
        pupilEllipse = [];
    end
    %for displaying completion
    if mod(iFrame,1000)==0
        disp([num2str(iFrame) ' / ' num2str(size(video,3)) ' frames']);
    end
    %define image frame and threshold
    im = video(:,:,iFrame);
    im_threshed = (im<threshDark).*maskEye;
    %rough estimate of center and size of pupil based of number of pixels
    [row1, col1] = find(im_threshed);
    %check to see if eye is open
    if length(row1)>threshNdarkPix
        pupilCenter_estimate = [mean(row1),mean(col1)];  %center
        pupilRadius_estimate = sqrt(length(row1)/pi);   %radius (from A=pi*r^2)
        %delete pixels which are too far away to be pupil
        for i=1:length(row1)
            %distance from estimated center
            dXY = pupilCenter_estimate - [row1(i), col1(i)];
            d = sqrt(sum(dXY.*dXY));
            %delete if too far
            if d>pupilRadius_estimate*radiusFudge
                im_threshed(row1(i),col1(i)) = false;
            end
        end
        %find edges in image with canny image detection
        %[~, threshold] = edge(im, 'canny');
        %im_edges = edge(im,'canny', threshold * edgeFudge).*maskEye;
        im_edges = edge(im,'canny').*maskEye;
        %delete edge pixels which are too far from center
        [row2, col2] = find(im_edges);
        for i=1:length(row2)
            %distance from estimated center
            dXY = pupilCenter_estimate - [row2(i), col2(i)];
            d = sqrt(sum(dXY.*dXY));
            %delete if too far
            if d>pupilRadius_estimate*radiusFudge
                im_edges(row2(i),col2(i)) = false;
            else
                %delete edge pixels which are too far from dark pixels
                [~,d] = dsearchn([row1,col1],[row2(i),col2(i)]);
                if d>threshPupilDist
                    im_edges(row2(i),col2(i)) = false;
                end
            end
        end
        %delete edge pixels which are not part of large group
        im_edges2 = bwlabel(im_edges); %label ROIs by unique number
        %number of ROIs in edges matrix
        nRois = max(im_edges2(:));
        %clear edge matrix and add back in if ROI is big
        im_edges = false(size(im_edges));
        for iROI = 1:nRois
            roiSize = sum(sum(im_edges2==iROI));
            if roiSize >= minEdgeSize
                im_edges(im_edges2==iROI) = true;
            end
        end
        %fit ellipse
        [row3,col3]=find(im_edges);
        pupilEllipse = fit_ellipse(row3,col3);
        %if ellipse found, use short axis as pupil diameter
        if ~isempty(pupilEllipse.long_axis) && pupilEllipse.long_axis~=0
            %record pupil diameter
            pupilD(iFrame) = pupilEllipse.long_axis;
            %record pupil position
            pupilXY(1,iFrame) = pupilEllipse.Y0_in;
            pupilXY(2,iFrame) = pupilEllipse.X0_in;
        end
    end
    
    %plot if you want to plot
    if plot_each_frame == true
        %plot pic with fitted ellipse
        pic = video(:,:,iFrame);
        for iii = 1:length(row3)
            pic(row3(iii),col3(iii)) = 255;
        end
        imshow(pic);
        %set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
        %pupil outline
        if ~isempty(pupilEllipse)
            ellipse(pupilEllipse.b,pupilEllipse.a,pupilEllipse.phi,pupilEllipse.Y0_in,pupilEllipse.X0_in,'r');
            %pupil center
            hold on
            scatter(pupilEllipse.Y0_in,pupilEllipse.X0_in,'r');
            hold off
        end
        pause;
        %F = getframe;
        %videoAnnotated(:,:,:,iFrame) = F.cdata;
    end
end
toc

close all

%% Delete pupil diameter measurements which are impossible

%any pupil diameters changing faster than this will be deleted
threshPupildRdT = 1.2;

pupilD_diff = abs(diff(pupilD));
[N,edges] = histcounts(pupilD_diff);
semilogy(edges(2:end),N)
hold on
line([threshPupildRdT threshPupildRdT], [1 max(N)],'Linewidth',2,'Color','r');
hold off
title('Changes in pupil diameter histogram')
ylabel('Count')
xlabel('Pixels/frame')

%% interpolate missing values

indDelete = find(pupilD_diff>threshPupildRdT)+1;

% figure
% plot(pupilD,'Linewidth',2)
% hold on
% scatter(indDelete,nanmean(pupilD)*ones(1,length(indDelete)),'filled')
% hold off

%delete values outside of dRdT range in pupilD and pupilXY
pupilD_inter = pupilD;
pupilXY_inter = pupilXY;
pupilD_inter(indDelete) = nan;
pupilXY_inter(:,indDelete) = nan;

pupilD_inter = inpaintn(pupilD_inter);
pupilXY_inter(1,:) = inpaintn(pupilXY(1,:));
pupilXY_inter(2,:) = inpaintn(pupilXY(2,:));

% figure
% plot(pupilD,'Linewidth',2)
% hold on
% plot(pupilD_inter,'Linewidth',2)
% scatter(indDelete,nanmean(pupilD)*ones(1,length(indDelete)),'filled','k')
% hold off

% calculate mean pupil position and delta position over time
pupilXY_mean = nanmean(pupilXY_inter,2);

%delta position over time
pupilXYDel = nan(1,length(pupilXY_inter));
for iFrame=1:length(pupilXY_inter)
    if ~isnan(pupilXY_inter(1,iFrame)) %avoid using values we deleted manually
        ydel = pupilXY_mean(1) - pupilXY_inter(1,iFrame);
        xdel = pupilXY_mean(2) - pupilXY_inter(2,iFrame);
        pupilXYDel(iFrame) = sqrt(ydel^2 + xdel^2);
    end
end

% imshow(mean(video,3),[0 255])
% hold on
% scatter(pupilXY_mean(1),pupilXY_mean(2));
% hold off
% 
% figure
% plot(pupilD_inter,'Linewidth',2);
% hold on
% plot(pupilXYDel+40,'Linewidth',2);
% hold off
% title('Pupil Diameter and Eye Position');
% xlabel('Time (1/10 sec)');
% ylabel('Pupil Diamter (pixels)');
% legend('Pupil diameter','Eye position');

% Low-pass filter pupil trace
%frequency cuttoff for low-pass filter
freqCutoff = 3;
sampleInterval = mode(diff(videoTS));
[b,a]=besself(4,2*pi*freqCutoff);
[bd,ad] = bilinear(b,a,1/sampleInterval);
pupilD_smooth = filtfilt(bd,ad,pupilD_inter-mean(pupilD_inter));
pupilD_smooth = pupilD_smooth + mean(pupilD_inter);

figure
plot(pupilD,'Linewidth',2)
hold on
plot(pupilD_smooth,'r','Linewidth',2)
title('Original and filtered pupil diameter measurements')
legend('Original', 'LP filtered');
hold off

correlogram = xcorr(pupilD_smooth,pupilD_inter,50,'coeff');
%plot(-50:1:50,correlogram)

%%
save('pupilD','pupilD')
save('pupilXY','pupilXY')
save('pupilXY_inter','pupilXY_inter')
save('pupilXY_mean','pupilXY_mean')
save('pupilXYDel','pupilXYDel')
save('pupilD_smooth','pupilD_smooth')
save('rectEye','rectEye')

%% Plot picture of eye and overlapping circle (optional)

figure
for iFrame = 1000:size(video,3)
    %plot video frame
    imshow(video(:,:,iFrame));
    %pupil outline
    ellipse(pupilD_smooth(iFrame)/2,pupilD_smooth(iFrame)/2,0,pupilXY_inter(1,iFrame),pupilXY_inter(2,iFrame),'r');
    %pupil center
    hold on
    scatter(pupilXY_inter(1,iFrame),pupilXY_inter(2,iFrame),'r');
    hold off
    pause;
end

%% Extra stuff for extra special people

figure
plot(wheelTime,smooth(wheelVelocity,11),'Linewidth',2)
xlabel('Time (seconds)');
ylabel('Wheel Speed (cm/s) and Pupil Diameter');
title('Mouse run speed');
hold on
%plot session, Task3
dotsize = 400;
yoffset = -5;
scatter(tStimAud,yoffset*ones(1,length(tStimAud)),dotsize,'m','filled');
scatter(tStimVis(logical(bIntHighVis)),yoffset*ones(1,sum(bIntHighVis)),dotsize,'r','filled','^');
scatter(tStimVis(logical(~bIntHighVis)),yoffset*ones(1,sum(~bIntHighVis)),dotsize,'r','filled','v');
%scatter(tStimVis,yoffset*ones(1,length(tStimVis)),dotsize,'r','filled');
scatter(tLicks,yoffset*ones(1,length(tLicks)),dotsize*.1,'k');
scatter(tResponses,yoffset*ones(1,length(tResponses)),dotsize*.7,'b');
%scatter(tStimCatch,yoffset*ones(1,length(tStimCatch)),dotsize,'b','filled');
%scatter(tStim,yoffset*bResponse(1:end))
plot(videoTS,pupilD_smooth*5-120,'r','Linewidth',2);
hold off

%% Cross correlation of pupil and running

maxlag = round(20 / .1);

%lp filter run speed
freqCutoff = 3;
sampleInterval = median(diff(wheelTime));
[b,a]=besself(4,2*pi*freqCutoff);
[bd,ad] = bilinear(b,a,1/sampleInterval);
%filter wheel velocity
wheelVelocitySmooth = filtfilt(bd,ad,wheelVelocity-mean(wheelVelocity))+mean(wheelVelocity);


%find common time matrix for wheel speed and pupil video
%but use deltaT of wheel speed
tStart = 4; %max(wheelTime(1),videoTS(1));
tEnd = 2270; %min(wheelTime(end),videoTS(end));
tCommon = wheelTime(find(wheelTime>=tStart,1,'first'):find(wheelTime<=tEnd,1,'last'));
%interpolate wheel data at common time
clipWheel = interp1(wheelTime(1:end-1),diff(wheelVelocitySmooth),tCommon);
%interpolate pupil data at common time
clipPupil = interp1(videoTS(1:end-1),diff(pupilD_smooth),tCommon);

%cross-correlate
[wheelPupilXcorr,lags] = xcorr(clipWheel-mean(clipWheel),clipPupil-mean(clipPupil),maxlag,'coeff');
tXcorr = lags*.1;

tMax = tXcorr(find(wheelPupilXcorr==max(wheelPupilXcorr)))

figure
plot(tXcorr,wheelPupilXcorr,'Linewidth',2)
hold on
line([0 0],[min(wheelPupilXcorr)-.05 max(wheelPupilXcorr)+.05],'Color','k');
title('dWheel/dt - dPupil/dt xcorr')
ylabel('Correlation')
xlabel('Lag (sec)')
ylim([-.12 .35])


