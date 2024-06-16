function [ewtc,mfb,Bw,Bt] = EWTC4Def(image)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % function [ewtc,mfb,Bw,Bt] = EWTC4Def(image)
    %
    % This function performs the Empirical Curvelet Transform for fringe 
    % separatio of transparent element.  
    %
    % Input:
    % image: input image (MUST BE SQUARE)
    % params: structure containing the following parameters:
    %
    % Output:
    % ewtc: cell containing each filtered output subband (ewtc{1} is the 
    %   lowpass subband the next ewtc{s}{t} are the bandpass filtered 
    %   images corresponds to the scales and t to the direction)
    % mfb: cell containing the set of empirical filters in the Fourier 
    %   domain (the indexation is the same as ewtc above)
    % Bw: list of the detected scale boundaries
    % Bt: list of the detected angle boundaries
    %
    % Author: Peide Yang / 202406 / Version 1.0
    %
    % Institution: Shanghai Engineering Research Center 
    %              of Ultra-Precision Optical Manufacturing, 
    %              School of Information Science and Technology, 
    %              Fudan University
    %
    % Algorithm is based on Empirical Wavelet Transforms code by Jerome 
    % Gilles and (Pseudo) Polar transform code by Michael Elad: 
    % https://ww2.mathworks.cn/matlabcentral/fileexchange/
    % 42141-empirical-wavelet-transforms
    %     
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
    %      [5] Averbuch, Amir, Ronald R. Coifman, David L. Donoho, Michael 
    %          Elad and Moshe Israeli. "Fast and accurate Polar Fourier 
    %          transform." Applied and Computational Harmonic Analysis 21 
    %          (2006): 145-167.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply a window function to the image and make sure the image is square
[H, W]=size(image);
originalImage = image;
Window = hamming(H) * hamming(W)';
image = image .* Window;
assert(H == W, 'Error: the image is not square');

%% Pseudo Polar FFT
image_PseudoFFT = PPFFT(image);

%% Defind the first scale (the low-pass filter)
indexppfft= IndexPPFFT(fftshift(abs(image_PseudoFFT')),'mean'); 
LL=round(length(indexppfft)/2);

boundariesW = 4; % it has to be an integer
Bw1 = boundariesW * pi / LL;

%% Find the anglar boundaries
indexppfft = IndexPPFFT(abs(image_PseudoFFT(floor(size(image_PseudoFFT,1)/2)+boundariesW(1):end,:)),'meanmax');
indexppfft = mapminmax(indexppfft,0,6000);

% Detect the angular boundaries use FTC_Seg method and find the top two maximum
boundaries = FTC_Seg([indexppfft,indexppfft,indexppfft],0);
boundaries = boundaries( ...
    (boundaries >= length(indexppfft)) & (boundaries <= 2 * length(indexppfft)));
boundaries = boundaries - length(indexppfft);
boundaries(boundaries == 0) = 1;
MaximumInTheInterval = zeros(1,length(boundaries)); 
for i =1:length(boundaries)-1
    [MaximumInTheInterval(i)] = max(indexppfft(boundaries(i):boundaries(i+1)));
    [MaximumInTheInterval(i+1)] = max([indexppfft(1:boundaries(1)),indexppfft(boundaries(end):end)]);
end
max_pos = zeros(1,length(MaximumInTheInterval));
for i = 1:length(MaximumInTheInterval)
    max_pos(i) = find(indexppfft == MaximumInTheInterval(i));
end
max_pos = sort(max_pos,2,"ascend");
MaximumInTheInterval = indexppfft(max_pos);
top_two_local_max_index = find_top_two_local_max_index_1D(MaximumInTheInterval); 
top_two_local_max_pos = max_pos(top_two_local_max_index);
corase_select_boundaries = findClosestSmallerIndices(boundaries, top_two_local_max_pos);
corase_select_boundaries_pos = boundaries(corase_select_boundaries);
corase_select_boundaries = unique([corase_select_boundaries,corase_select_boundaries+1]); 
corase_select_boundaries(corase_select_boundaries > length(boundaries)) = 1; 
corase_select_boundaries = sort(corase_select_boundaries,'ascend');
boundaries = unique(boundaries(corase_select_boundaries));
corase_select_boundaries_pos(indexppfft(corase_select_boundaries_pos) > 2000) = [];
boundaries(indexppfft(boundaries) > 2000) = [];



% figure;
% hold on;
% title('anglar boundaries');
% plot(indexppfft); 
% for k = 1:length(boundaries)
%     plot([boundaries(k) boundaries(k)], [min(indexppfft) max(indexppfft)], '--r'); 
% end
% hold off;

%% Compute the meanmax spectrum with respect to the angles to find the radial boundaries
% we detect the scales per angular sector , only detect the sectors which
% have the maximum
Bt = (boundaries-1)*pi/length(indexppfft)-3*pi/4;
Bw=cell(length(boundaries)+1,1);
Bw{1}=Bw1;
select_angle_num = 0;
for t=1:length(boundaries)-1
    if any(boundaries(t)==corase_select_boundaries_pos)
    % spectrum on the given angular range
    indexppfft = IndexPPFFT(abs(image_PseudoFFT((boundariesW(1)+ ...
        floor(size(image_PseudoFFT,1)/2)):end,boundaries(t):boundaries(t+1)))','meanmax');
    indexppfft = mapminmax(indexppfft,0,10000);
    % Detect the boundaries in radial direction
    bounds = FTC_Seg(indexppfft,0);
    if length(bounds)>=2 && any(bounds>300)
        [max_of_the_last_mode,max_pos_of_the_last_mode] = max(indexppfft(bounds(end-1):bounds(end)));
        denoiseIndex = find(indexppfft(bounds(end-1)+max_pos_of_the_last_mode:end) < 0.01*max_of_the_last_mode ...
            , 1, 'first'); % Limit the radial boundaries
        bounds = sort([bounds,bounds(end-1)+max_pos_of_the_last_mode+denoiseIndex+10]);
    end
    bounds = abs(bounds);
    bounds(bounds > 300)=300;
    Bw{t+1} = (boundariesW(1) + bounds)*pi/LL;
%     figure;
%     hold on;
%     title('radial boundaries');
%     plot(indexppfft); 
%     for k = 1:length(bounds)
%         plot([bounds(k) bounds(k)], [min(indexppfft) max(indexppfft)], '--b'); 
%     end
%     hold off;
    select_angle_num = select_angle_num + 1;
    end
end

% last one 
if select_angle_num < 2
    indexppfft=sum(abs(image_PseudoFFT(( ...
    boundariesW(1)+floor(size(image_PseudoFFT,1)/2)):end, ...
    boundaries(end):end)),2)+sum(abs(image_PseudoFFT(( ...
    boundariesW(1)+floor(size(image_PseudoFFT,1)/2)):end, ...
    1:boundaries(1))),2);
    % Detect the boundaries
    indexppfft = mapminmax(indexppfft',0,10000);
    bounds = FTC_Seg(indexppfft,0);
    if length(bounds)>=2 
        [max_of_the_last_mode,max_pos_of_the_last_mode] = max(indexppfft(bounds(end-1):bounds(end)));
        denoiseIndex = find(indexppfft(bounds(end-1)+max_pos_of_the_last_mode:end) < 0.01*max_of_the_last_mode ...
            , 1, 'first');
        bounds = sort([bounds,bounds(end-1)+max_pos_of_the_last_mode+denoiseIndex+10]);
        if any(bounds>length(indexppfft))
            bounds(bounds>length(indexppfft)) = length(indexppfft);
        end
    end
%     figure;
%     hold on;
%     title('radial boundaries');
%     plot(indexppfft); 
%     for k = 1:length(bounds)
%         plot([bounds(k) bounds(k)], [min(indexppfft) max(indexppfft)], '--b');
%     end
%     hold off;

    Bw{end} = (boundariesW(1)+bounds)*pi/LL;
end

% Set an upper bound for the empty Bw
drop_idx = 0;
for i = 2:length(Bw) 
    if isempty(Bw{i})
        Bw{i} = 3;
        drop_idx = [drop_idx,i];
    else
        Bw{i} = unique(Bw{i});
    end
end

% Build the filter bank
mfb = EWT2D_Curvelet_FilterBank(Bw,Bt,W,H,3);

% Filter the original image to extract each subband
ff=fft2(originalImage); 

% Extract the low frequencies first
ewtc = cell(length(mfb),1);
ewtc{1} = ifft2(conj(mfb{1}).*ff);

% Extract the other frequencies
n = 1;
for s=2:length(mfb)
    if s ~= drop_idx
        for t=1:length(mfb{s})
            ewtc{n} = real(ifft2(conj(mfb{s}{t}).*ff));
            n = n + 1;
        end
    end
end  
end