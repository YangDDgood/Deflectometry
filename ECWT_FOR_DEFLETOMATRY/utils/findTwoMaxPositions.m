function [max_positions] = findTwoMaxPositions(arr)
    [~, max_index1] = max(arr);
    arr(max_index1) = -inf; 
    [~, max_index2] = max(arr);
    max_positions = [max_index1, max_index2];
end