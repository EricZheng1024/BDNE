function [L, P] = sub_rects_3(point, point_2, L, P, i_point_cut, is_obj)
% 只切必要的矩形
% 矩形是否与point和point_2规定的矩形有相交
% P: 矩形顶点的目标值为解的目标值的逻辑索引
% age: 矩形年龄
% is_obj: point_cut是否为目标向量

    switch i_point_cut
        case 1
            point_cut = point;
        case 2
            point_cut = point_2;
        otherwise
            error('Incorrect index of points for cutting.');
    end
    for s = 1 : length(L)
        Rs = L{1};
        C_s_prime = Rs(:,1)<point_cut & point_cut<Rs(:,2);
        if any(C_s_prime) && ...
                all( (point<=Rs(:,1) & Rs(:,1)<=point_2) | (point<=Rs(:,2) & Rs(:,2)<=point_2) ...
                | (Rs(:,1)<=point & point<=Rs(:,2)) | (Rs(:,1)<=point_2 & point_2<=Rs(:,2)) )  % 必须等号
            C_s_prime = find(C_s_prime);
            L_prime = {Rs};
            P_prime = P(1);
            for j = C_s_prime'
                L_pprime = cell(1, 2*length(L_prime));
                P_pprime = cell(1, 2*length(P_prime));
                for t = 1 : length(L_prime)
                    %
                    R1 = L_prime{t};
                    R1(j,2) = point_cut(j);
                    L_pprime{2*t-1} = R1;
                    % L_pprime{t} = R1;
                    %
                    Q1 = P_prime{t};
                    Q1(j,2) = is_obj;
                    P_pprime{2*t-1} = Q1;
                    % P_pprime{t} = Q1;
                    %
                    R2 = L_prime{t};
                    R2(j,1) = point_cut(j);
                    L_pprime{2*t} = R2;
                    % L_pprime{t+length(L_prime)} = R2;
                    %
                    Q2 = P_prime{t};
                    Q2(j,1) = is_obj;
                    P_pprime{2*t} = Q2;
                    % P_pprime{t+length(L_prime)} = Q2;
                end
                L_prime = L_pprime;
                P_prime = P_pprime;
            end
            L = [L(2:end) L_prime]; P = [P(2:end) P_pprime];
        else
            L = [L(2:end) L(1)]; P = [P(2:end) P(1)];
        end
    end
end
