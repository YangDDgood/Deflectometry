%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a demo to perform the Empirical Curvelet Transform for fringe 
% separation of transparent elements.
%
% Author: Peide Yang / 202406 / Version 1.0
%
% Institution: Shanghai Engineering Research Center 
%              of Ultra-Precision Optical Manufacturing, 
%              School of Information Science and Technology, 
%              Fudan University
%
% Algorithm is based on Empirical Wavelet Transforms code by Jerome Gilles
% and (Pseudo) Polar transform code by Michael Elad: 
% https://ww2.mathworks.cn/matlabcentral/fileexchange/
% 42141-empirical-wavelet-transforms
% Ref. [1] Recognition and separation of fringe patterns in 
%          deflectometric measurement of transparent elements based on 
%          empirical curvelet transform 
%      [2] J.Gilles, "Empirical wavelet transform", IEEE Trans. 
%          Signal Processing, 2013.
%      [3] J.Gilles, G.Tran, S.Osher "2D Empirical transforms. 
%          Wavelets, Ridgelets and Curvelets Revisited", SIAM Journal 
%          on Imaging Sciences, Vol.7, No.1, 157--186, 2014.
%      [4] J. Delon, A. Desolneux, J-L. Lisani and A-B. Petro, A non
%          parametric approach for histogram segmentation, IEEE 
%          Transactions on Image Processing, vol.16, no 1, pp.253-261, 
%          Jan. 2007. 
%      [5] Averbuch, Amir, Ronald R. Coifman, David L. Donoho, Michael Elad
%          and Moshe Israeli. "Fast and accurate Polar Fourier transform." 
%          Applied and Computational Harmonic Analysis 21 (2006): 145-167.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clf;close all;clc; 
addpath '.\utils'
addpath '.\Empirical Wavelet Transforms'
addpath '.\experiment'
addpath '.\simulation'

img = mat2gray(im2double(imread('./experiment/2.bmp')));
[H,W] = size(img);
padsize = 24;
img = padarray(img, [padsize,padsize], "replicate",'both');
Window = hamming(size(img,1)) * hamming(size(img,2))';
s_image = img .* Window;
ff=log(1+fftshift(abs(fft2(s_image))));
figure;
subplot(121);imagesc(s_image);
title('the image used for filter bank construction');colorbar;colormap gray;
subplot(122);imagesc(ff);
title('Fourier transform of the windowed image');colorbar; colormap gray ;
%% perform ECWT for deflectometry
[ewtc,mfb,Bw,Bt] = EWTC4Def(img);
for i = 1: length(ewtc)
    ewtc{i} = ewtc{i}(padsize+1:end-padsize, padsize+1:end-padsize);
end

Show_Curvelets_boundaries(img,Bw,Bt,3)
figure;
cn = ceil(sqrt(length(ewtc)+1));
cn2 = ceil((length(ewtc)+1)/cn);
subplot(cn2,cn,1); imagesc(img);colormap gray;colorbar;title('Input image');
for tt = 2:length(ewtc)+1
    subplot(cn2,cn,tt); imagesc(ewtc{tt-1}); colormap gray; colorbar;
    title(['all the Modes ',num2str(tt-1),'/',num2str(length(ewtc))]); 
end
%% Calculate weighted permutation entropy to remove noise modes
ent_c = zeros(1,length(ewtc));
ent_r = zeros(1,length(ewtc));
ent_min = zeros(1,length(ewtc));
for i = 1:length(ewtc)
    ent_r(i) = getPermutationEntropy(ewtc{i}, 5, 1);
    ent_c(i) = getPermutationEntropy(ewtc{i}', 5, 1);
    ent_min(i) = min(ent_r(i), ent_c(i));
end
% Set a threshold to delete noise
ewtc = ewtc(ent_min > 1e-4 & ent_min < 0.5); 

%% Get variance and initial fringes
variance =zeros(1,length(ewtc));
for i = 1:length(ewtc)
    variance(i) = var(ewtc{i},1,'all');
end

% Further determine meaningful modes
waterLine = watershed(-variance); 
[max_positions] = findTwoMaxPositions(variance);
region1 = (waterLine == waterLine(max_positions(1)));
region2 = (waterLine == waterLine(max_positions(2)));
region1 = imerode(region1, strel('line', 3, 0));
region2 = imerode(region2, strel('line', 3, 0));
region1(max_positions(2)) = 0;
region2(max_positions(1)) = 0;
ewtcImage1 = ewtc(region1); 
ewtcImage2 = ewtc(region2); 
initialModes = ewtc(max_positions);

figure; 
cn = ceil(sqrt(length(ewtcImage1)+length(ewtcImage2)+1));
cn2 = ceil((length(ewtcImage1)+length(ewtcImage2)+1)/cn);
subplot(cn2,cn,1); imagesc(img); colormap gray; colorbar; 
title('Input image');
for tt = 2 : length(ewtcImage1) + 1
    subplot(cn2, cn, tt); imagesc(ewtcImage1{tt-1});colormap gray;colorbar;
    title(['meaningfulModes ',num2str(tt-1),'/', ...
    num2str(length(ewtcImage1)+length(ewtcImage2))]); 
end
for tt = length(ewtcImage1) + 2 : length(ewtcImage1)+length(ewtcImage2)+1
    subplot(cn2, cn, tt); imagesc(ewtcImage2{tt-length(ewtcImage1)-1});
    colormap gray;colorbar;title(['meaningfulModes ',num2str(tt-1), '/',...
    num2str(length(ewtcImage1)+length(ewtcImage2))]); 
end

%% combination
[outputImage1, combinationIdx1] = getModesConbination(initialModes{1}, ewtcImage1);
[outputImage2, combinationIdx2] = getModesConbination(initialModes{2}, ewtcImage2);
%% Map the images to 0-1 (optional)
if min(outputImage1(:)) < 0
    outputImage1 = outputImage1 + abs(min(outputImage1(:)));
end
if min(outputImage2(:)) < 0
    outputImage2 = outputImage2 + abs(min(outputImage2(:)));
end
if max(outputImage1(:)) > 1
    outputImage1 = outputImage1 / max(outputImage1(:));
end
if max(outputImage2(:)) > 1
    outputImage2 = outputImage2 / max(outputImage2(:));
end
outputImage1 = imadjust(outputImage1, stretchlim(outputImage1, [0 1]), [0 1]);
outputImage2 = imadjust(outputImage2, stretchlim(outputImage2, [0 1]), [0 1]);
% Show Fig
figure;imagesc(outputImage1);colormap gray;title("outputImage1");
figure;imagesc(outputImage2);colormap gray;title("outputImage2");
% Save
imwrite(outputImage1, fullfile('.\', ['ouput1', '.bmp']));
imwrite(outputImage2, fullfile('.\', ['ouput2', '.bmp']));