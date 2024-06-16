function WPE = getPermutationEntropy(x, m, tau)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the weighted permutation entropy of a time series.
    %
    % WPE = getPermutationEntropy(x, m, tau) 
    %     
    % Inputs:
    % m: Embedding Dimension
    % tau: Time Delay
    % Output:
    % WPE: weighted permutation entropy
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
    %      [2] Bilal Fadlallah, Badong Chen, Andreas Keil, et al. 
    %          Weighted-permutation entropy: A complexity measure for 
    %          time series incorporating amplitude information[J]. 
    %          Physical review. E, Statistical, nonlinear, 
    %          and soft matter physics, 2013, 87(2): 022911. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x = squeeze(x);
    x = x(:);
    x = single(x);% for efficient
    len = length(x);
    cols = len - (m-1) * tau;
    dataMat = zeros(m, cols, "single");
    if m < cols
        for i = 1:m
            dataMat(i,:) = x((i-1)*tau+1:(i-1)*tau+cols);
        end
    else
        for i = 1:cols
            dataMat(:,i) = x(i:tau:i+m*tau-1);
        end
    end
    % Compute weight
    w = var(dataMat,1); 
    w = w/sum(w);
    [~,s] = sort(dataMat); 
    s = s';
    S = unique(s,'rows'); 
    p = zeros(size(S,1),1);
    for i = 1:size(S,1)
        ind = sum(abs(s-S(i,:)),2)==0;
        p(i) = sum(w(ind));
    end
    WPE = -sum(p.*log(p))/log(factorial(m)); 
end