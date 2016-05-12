clear;clc;
ImgNo = 07;
% The only input parameter is the image number.
% However, the parameters depend on the image we choose.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Explanation on these parameters:
% ScaleH, ScaleL: Crop the original image to exclude nonsense pixels. In
%                 other words, we just focus on the road.
% PixelGap      : Reduce the resolution of this image to exclude the 
%                 outlines of small size.
% Thresh_Edge   : Used in the edge function. It's like a filter. Smaller
%                 values result in edges of smaller size.
% RhoRes_Hough  : Used in the hough function. Smaller values lead more
%                 careful hough transformation.
% Theta_Hough   : Used in the hough function. Smaller ranges would avoid
%                 parallel lines to be found.
% Thresh_Peak   : Used in the houghpeaks function. Smaller values cause the
%               : function to get more lines.
% PeakN         : Used in the houghpeaks function. Smaller values give less
%                 lines.
% FillGap       : Used in the houghlines function. Combine lines accoording
%                 to this value. Larger values cause lines to be longer.
% MinLen        : Used in the houghlines function. The length of lines
%                 found must be larger than this value. Smaller values give
%                 more lines.
%NHoodSize_Scale: Used in the houghpeaks function. It seems that smaller
%                 values give more lines.
switch ImgNo
    case {01, 02}
        ScaleH = 0.75;
        ScaleL = 0.65;
        PixelGap = 1;
        Thresh_Edge = 0.3;
        RhoRes_Hough = 0.3;
        Theta_Hough = -80:0.5:80;
        Thresh_Peak = 0.3;
        PeakN = 6;
        FillGap = 40;
        MinLen = 10;
        NHoodSize_Scale = 30;
%     case 02
%         ScaleH = 0.75;
%         ScaleL = 0.65;
%         PixelGap = 1;
%         Thresh_Edge = 0.3;
%         RhoRes_Hough = 0.3;
%         Theta_Hough = -80:0.5:80;
%         Thresh_Peak = 0.3;
%         PeakN = 6;
%         FillGap = 40;
%         MinLen = 10;
    case {03, 04}
        ScaleH = 1.0;
        ScaleL = 0.6;
        PixelGap = 3;
        Thresh_Edge = 0.3;
        RhoRes_Hough = 0.3;
        Theta_Hough = -80:0.5:80;
        Thresh_Peak = 0.3;
        PeakN = 6;
        FillGap = 10;
        MinLen = 30;
        NHoodSize_Scale = 30;
%     case 04
%         ScaleH = 1.0;
%         ScaleL = 0.6;
%         PixelGap = 3;
%         Thresh_Edge = 0.3;
%         RhoRes_Hough = 0.3;
%         Theta_Hough = -80:0.5:80;
%         Thresh_Peak = 0.3;
%         PeakN = 6;
%         FillGap = 10;
%         MinLen = 30;
    case 05 %
        % The centerline should be detected by other methods.
        ScaleH = 1.0;
        ScaleL = 0.6;
        PixelGap = 2;
        Thresh_Edge = 0.5;
        RhoRes_Hough = 0.1;
        Theta_Hough = -85:0.5:85;
        Thresh_Peak = 0.1;
        PeakN = 6;
        FillGap = 5;
        MinLen = 5;
        NHoodSize_Scale = 30;
    case 06
        ScaleH = 1.0;
        ScaleL = 0.6;
        PixelGap = 2;
        Thresh_Edge = 0.5;
        RhoRes_Hough = 0.5;
        Theta_Hough = -85:0.5:85;
        Thresh_Peak = 0.5;
        PeakN = 6;
        FillGap = 10;
        MinLen = 10;
        NHoodSize_Scale = 30;
    case 07
        ScaleH = 0.95;
        ScaleL = 0.7;
        PixelGap = 2;
        Thresh_Edge = 0.4;
        RhoRes_Hough = 0.5;
        Theta_Hough = -65:0.5:65;
        Thresh_Peak = 0.3;
        PeakN = 6;
        FillGap = 40;
        MinLen = 10;
        NHoodSize_Scale = 30;
    case 08
        ScaleH = 1.0;
        ScaleL = 0.6;
        PixelGap = 2;
        Thresh_Edge = 0.5;
        RhoRes_Hough = 0.5;
        Theta_Hough = -85:0.5:85;
        Thresh_Peak = 0.5;
        PeakN = 6;
        FillGap = 40;
        MinLen = 30;
        NHoodSize_Scale = 30;
end
% Read in image
ImgPath = strcat('./img/', sprintf('%02d', ImgNo));
ImgPath = strcat(ImgPath, '.jpg');
Img_RGB = imread(ImgPath);
Img_Size = size(Img_RGB);
% RGB to Gray
Img_Gray = rgb2gray(Img_RGB);
% Crop the image
ImgTarget = Img_Gray(round(Img_Size(1)*ScaleL):round(Img_Size(1)*ScaleH), :, :);
% Filter pixels
I = ImgTarget(1:PixelGap:end, 1:PixelGap:end);
figure
PlotRowN = 2;
PlotColN = 1;
subplot(PlotRowN,PlotColN,1);
imshow(I);

BW = edge(I, 'canny', Thresh_Edge);
subplot(PlotRowN,PlotColN,2); hold on;
imshow(BW);
% Plot the location of the observer.
BW_Size = size(BW);
PtsCenterX = round(BW_Size(2)/2);
PtsCenterY = BW_Size(1);
plot(PtsCenterX, PtsCenterY, 'o','LineWidth',2,'Color','blue');

% Hough transform
[H, theta, rho] = hough(BW, ...
    'RhoResolution', RhoRes_Hough, 'Theta', Theta_Hough);
% figure, imshow(imadjust(mat2gray(H)),[],'XData',theta,'YData',rho,...
%         'InitialMagnification','fit');
% xlabel('\theta (degrees)'), ylabel('\rho');
% axis on, axis normal, hold on;
% colormap(hot)
% Find the votes peaks in H, which represents the target lines.
NHoodSize = ceil(size(H)/NHoodSize_Scale);
NHoodSize = NHoodSize + ~mod(NHoodSize, 2);
P = houghpeaks(H, PeakN, ...
    'Threshold', ceil(Thresh_Peak*max(H(:))), ...
    'NHoodSize', NHoodSize);
% P = houghpeaks(H, PeakN, ...
%     'Threshold', ceil(Thresh_Peak*max(H(:))) );
% Find the found lines in the real gray image
lines = houghlines(BW,theta,rho,P, ...
    'FillGap', FillGap, 'MinLength', MinLen);
% Plot the lines
EndPts = zeros(length(lines), 2, 2);
lines_length = zeros(length(lines), 1);
distances = zeros(length(lines), 1);
angles_ac = zeros(length(lines), 1);
angles_hor = zeros(length(lines), 1);
heights = zeros(length(lines), 1);

for k = 1:length(lines)
    EndPts(k,:,:) = [lines(k).point1; lines(k).point2];
    plot(EndPts(k,:,1),EndPts(k,:,2),'LineWidth',2,'Color','green');
    % Plot beginnings and ends of lines
    plot(EndPts(k,1,1),EndPts(k,1,2),'x','LineWidth',2,'Color','yellow');
    plot(EndPts(k,2,1),EndPts(k,2,2),'x','LineWidth',2,'Color','red');
    c = [EndPts(k,1,1)-EndPts(k,2,1), EndPts(k,1,2)-EndPts(k,2,2)];
    c_mag = sqrt((c(1))^2 + (c(2))^2);
    a = [EndPts(k,1,1)-PtsCenterX, EndPts(k,1,2)-PtsCenterY];
    a_mag = sqrt((a(1))^2 + (a(2))^2);
    lines_length(k) = c_mag;
    distances(k) = a_mag;
    angles_ac(k) = acos(dot(a, c) / (a_mag * c_mag));
    angles_hor(k) = atan(c(2)/c(1));
    heights(k) = a_mag * sin(angles_ac(k));
%     disp(k)
%     disp(lines_length(k))
%     disp(angles_hor(k) / pi * 180)
%     disp(heights(k))
%     pause
end

% Divide the whole series of lines into 2 classes, i.e. Left and Right, in
% terms of the angle between the line and the horizon.
anglesL_Ind = find(angles_hor*180/pi <= 0);
anglesR_Ind = find(angles_hor*180/pi > 0);
[LineL_height, LineL_Ind] = min(heights(anglesL_Ind));
[LineR_height, LineR_Ind] = min(heights(anglesR_Ind));
k = anglesL_Ind(LineL_Ind);
plot(EndPts(k,:,1),EndPts(k,:,2),'LineWidth',2,'Color','blue');
disp('Angle of the left line')
disp(angles_hor(k)*180/pi);
k = anglesR_Ind(LineR_Ind);
plot(EndPts(k,:,1),EndPts(k,:,2),'LineWidth',2,'Color','blue');
disp('Angle of the right line')
disp(angles_hor(k)*180/pi);

ImgNewPath = ImgPath(1:end-4);
saveas(gcf, ImgNewPath, 'png');
