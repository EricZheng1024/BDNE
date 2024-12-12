function score = DistNadir(Population,optimum,p)
% <min>
% Lp-distance of estimated nadir point to nadir point

%------------------------------- Reference --------------------------------
% 
%--------------------------------------------------------------------------

    PopObj = Population.best.objs;
    if size(PopObj,2) ~= size(optimum,2)
        score = nan;
    else
        switch nargin
            case 2
                p = 2;
            case 3
                % do noting
            otherwise
                error('The number of input parameters is wrong.')
        end
        pop_max = max(PopObj, [], 1);
        opt_min = min(optimum, [], 1);
        opt_max = max(optimum, [], 1);
        score = norm((pop_max-opt_max)./(opt_max-opt_min), p);
    end
end