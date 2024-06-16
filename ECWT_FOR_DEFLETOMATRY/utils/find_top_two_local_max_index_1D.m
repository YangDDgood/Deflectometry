
function output = find_top_two_local_max_index_1D(vector)
    max1 = -Inf;
    max2 = -Inf;
    index1 = -1;
    index2 = -1;

    for i = 1:length(vector)
        num = vector(i);
        if num > max1
            max2 = max1;
            index2 = index1;
            max1 = num;
            index1 = i;
        elseif num > max2 && num < max1
            max2 = num;
            index2 = i;
        end
    end
output = [index1, index2];
end