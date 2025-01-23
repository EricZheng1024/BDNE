function [L, is_remove] = remove_rects_2(L, l, u, tol)
% 返回is_remove

    is_remove = false(length(L),1);
    for i = 1 : length(L)
        % if all(L{i}(:,1)>=l) && all(L{i}(:,2)<=u)  % Theoretically, it is equivalent to any(L{i}(:,1)>=l) && any(L{i}(:,2)<=u) since rectangles are non-overlapping. 
        if all(L{i}(:,1)-l>=-tol) && all(L{i}(:,2)-u<=tol)
            is_remove(i) = true;
        end
    end
    L(is_remove) = [];
end