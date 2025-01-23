function V = cal_rect_vols(L)
    V = zeros(length(L),1);
    for i = 1 : length(L)
        tmp = diff(L{i},[],2);
        % if any(tmp<0)  % For debugging
        %     error('The lower vertex of the rectangle does not dominate the upper vertex.')
        % end
        V(i) = prod(tmp);
    end
    % if any(V<1e-15)  % For debugging
    %     error('The volume of rectangle should be large than zero.')
    % end
end