function [output_image, combinationIdx] = getModesConbination(initialMode, modes, th)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this function is used to conbine the possible over-segmentation modes
% input: 
% initialMode : the initial fringe mode
% modes : The possible missing part of the initial mode
% th :threshold (if necessary)
% output_image : the output image
% 
% Author: Peide Yang / 202406 / Version 1.0
%
% Institution: Shanghai Engineering Research Center 
%              of Ultra-Precision Optical Manufacturing, 
%              School of Information Science and Technology, 
%              Fudan University
%
% ref. Recognition and separation of fringe patterns in deflectometric 
% measurement of transparent elements based on empirical curvelet transform 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if nargin == 2
        th = 0;
    end
    output_image = initialMode;
    combinationIdx = [];
    Variance =zeros(1,length(modes));
    for i =1:length(modes)
        Variance(i) = var(modes{i}, 1, 'all');
    end

    for i = 1:length(modes)
        entOldr = getPermutationEntropy(output_image, 4, 1);
        entOldc = getPermutationEntropy(output_image', 4, 1);
        entOld = min(entOldr, entOldc);
        
        entNewr = getPermutationEntropy( output_image + modes{i}, 4, 1);
        entNewc = getPermutationEntropy((output_image + modes{i})', 4, 1);
        entNew = min(entNewr, entNewc);
        if entOld > entNew && (Variance(i) > th)
            output_image = output_image + modes{i};
            combinationIdx = [combinationIdx, i];
        end
    end
end