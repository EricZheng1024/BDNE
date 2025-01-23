function [L, V] = sub_rects_2(point, point_2, L, V, point_cut)
% 只切必要的矩形
% 矩形是否与point和point_2规定的矩形有相交

    switch nargin
        case 4
            point_cut = point;
        case 5
            % do noting
        otherwise
            error('Incorrect number of input variables.');
    end
    for s = 1 : length(L)
        Rs = L{1};
        C_s_prime = Rs(:,1)<point_cut & point_cut<Rs(:,2);
        if any(C_s_prime) && ...
                all( (point<=Rs(:,1) & Rs(:,1)<=point_2) | (point<=Rs(:,2) & Rs(:,2)<=point_2) ...
                | (Rs(:,1)<=point & point<=Rs(:,2)) | (Rs(:,1)<=point_2 & point_2<=Rs(:,2)) )  % 必须等号
            C_s_prime = find(C_s_prime);
            L_prime = {Rs};
            for j = C_s_prime'
                L_pprime = cell(1, 2*length(L_prime));
                for t = 1 : length(L_prime)
                    R1 = L_prime{t};
                    R1(j,2) = point_cut(j);
                    L_pprime{2*t-1} = R1;
                    R2 = L_prime{t};
                    R2(j,1) = point_cut(j);
                    L_pprime{2*t} = R2;
                end
                L_prime = L_pprime;
            end
            L = [L(2:end) L_prime]; V = [V(2:end); cal_rect_vols(L_prime)];
        else
            L = [L(2:end) L(1)]; V = [V(2:end); V(1)];
        end
    end
end
