function index = IndexPPFFT(PPFFT, index_type)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculates the indicated value of PPFFT result
    % index_of_everyangle = IndexPPFFT(PPFFT, index_type)
    %     
    % Inputs:
    % PPFFT: the PPFFT result
    % index_type:  the methods to get index for angular and radial sectors 
    % segmentation
    % Output:
    % index: the 1D index for angular and radial sectors
    %     
    % Author: Peide Yang / 202406 / Version 1.0
    %
    % Institution: Shanghai Engineering Research Center 
    %              of Ultra-Precision Optical Manufacturing, 
    %              School of Information Science and Technology, 
    %              Fudan University
    %
    % Ref. [1] Recognition and separation of fringe patterns in 
    %          deflectometric measurement of transparent elements based on 
    %          empirical curvelet transform 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~, N] = size(PPFFT);
    if strcmp(index_type,'mean')
        index = mean(PPFFT,1);
    elseif strcmp(index_type,'max')
        index = max(PPFFT,1);
    elseif strcmp(index_type,'meanmax')
        index = zeros(1,N);
        for k = 1:N
        elements_in_everyangle = sort(PPFFT(:,k), 'descend');
        top1PercentPos = ceil(numel(elements_in_everyangle) * 0.1);
        % Take the elements before the first 0.5% position and 
        % calculate the mean of the first 0.5% elements
        index(1,k) = mean(elements_in_everyangle(1:top1PercentPos));
        end
    else 
        error('unknown index_type')
    end
end