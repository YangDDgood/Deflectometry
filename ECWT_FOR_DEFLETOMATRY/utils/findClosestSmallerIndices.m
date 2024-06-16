function closestIndices = findClosestSmallerIndices(b, a)
    closestIndices = zeros(size(a));  
    for i = 1:numel(a)
        smallerIndices = find(b < a(i));
        if ~isempty(smallerIndices)
            [~, idx] = max(smallerIndices);
            closestIndices(i) = smallerIndices(idx);
        else
            closestIndices(i) = length(b);
        end
    end
end