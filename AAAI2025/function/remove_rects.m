function [L, V] = remove_rects(L, V, l, u, tol)
    is_remove = false(length(L),1);
    for i = 1 : length(L)
        % if all(L{i}(:,1)>=l) && all(L{i}(:,2)<=u)  % Theoretically, it is equivalent to any(L{i}(:,1)>=l) && any(L{i}(:,2)<=u) since rectangles are non-overlapping. 
        if all(L{i}(:,1)-l>=-tol) && all(L{i}(:,2)-u<=tol)
            is_remove(i) = true;
        end
    end
    L(is_remove) = [];
    V(is_remove) = [];
end
