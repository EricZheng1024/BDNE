% clear, close all

function [y_N,Y_P,n_model_solved,runtime,runtime_solver,runtime_k] = BDNE_MOIQP(H,f,A,b,lb,ub,Aeq,beq,min_finterval,p, ...
    filename, ...
    runtime0_k,options)

%% Config
% general_config;
mu = 1e6;
tol = min_finterval / 2;  % 删除矩形的容差

n_model_solved = zeros(p, 1);
runtime_solver = 0;
runtime_k = zeros(p,2);
start_time = cputime;


%% Ideal objective vector
[y_I, n_model_solved, runtime_solver] = payoff_table_MOIQP(H,f,A,b,lb,ub,Aeq,beq,p, ...
    options, n_model_solved, runtime_solver);


%% Upper objective values
[y_U, n_model_solved, runtime_solver] = upper_bound_MOIQP(H,f,A,b,lb,ub,Aeq,beq,p, ...
    options, n_model_solved, runtime_solver);


%% Main procedure
runtime_pre = cputime-start_time;
Y_P = cell(1, p);
C = cell(1, p);
y_N = y_I;
% First guess of range of normalization  若z_r1与z_r2改变，注意更新相关参数  范围必须包裹整个PF，否则？
z_r1 = y_I;
z_r2 = y_U;
% Fixed parameters of the boundary subporblem
% gain_norm = 1e4;  % 得考虑归一化范围、alpha的设置、有效位数（15~17）
% gain_weight = 1e4;
gain_norm = 10^ceil(log10(max(y_U-y_I)));  % 如果目标函数范围是小数，则不适用
if mu == inf
    % alpha = 0;
    alpha_gain = 0;
    minus_alpha_gain = 1;
else
    % alpha = p / ((p-1) * mu + p);  % 太小会使得intlinprog出现精度问题，使得结果有问题，而且不会报错！
    tmp = 10^max(8-log10(gain_norm),0);
    tmp = 1;
    alpha_gain = p * tmp / ((p-1) * mu + p);
    minus_alpha_gain = tmp - alpha_gain;
end
f_bd = zeros(size(f,2)+1,1);
f_bd(1) = 1;
Aeq_bd = [zeros(size(Aeq,1),1)  Aeq];
beq_bd = beq;
lb_bd = [0; lb];
ub_bd = [inf; ub];
% tol = trans_cone(tol./(z_r2-z_r1), alpha);
% tol = norm_trans(tol, 0, z_r2-z_r1, alpha, gain_norm);
% tol = norm_trans(tol, 0, z_r2-z_r1, alpha_gain, minus_alpha_gain, gain_norm);
% tol_precision = 1e-5*gain_norm;

% options_fea = options;
% options_fea.MaxFeasiblePoints = 1;
% options_fea.MaxTime = 2;
% options_fea.NodeSelection = 'mininfeas';
% options_fea.CutGeneration = 'none';
% % options_fea.CutMaxIterations = 1;
% options_fea.Heuristics = 'none';
% % options_fea.HeuristicsMaxNodes = 1;
% options_fea.BranchRule = 'mostfractional';  % 由于目标函数为0，所以伪代价无法计算
% % options_fea.LPPreprocess = 'none';
% % options_fea.IntegerPreprocess = 'advanced';
% % options_fea.RootLPAlgorithm = "primal-simplex";
options_fea = options;
options_fea.solverOpts = {'limits/solutions',1; ...
                          'limits/time',2; ...
                          'numerics/feastol',1e-7; ...
                          'numerics/epsilon',1e-9; ...
                          'limits/gap',0; ...
                          'limits/absgap',0};

% For each objective
for k = 1 : p
    % Init
    start_time_k = cputime;
    C{k} = setdiff(1:p,k);  % 从小到大
    % z_r1_tilde = trans_cone(normalize_2(z_r1, z_r1, z_r2), alpha);
    % z_r2_tilde = trans_cone(normalize_2(z_r2, z_r1, z_r2), alpha);
    % z_r1_tilde = norm_trans(z_r1, z_r1, z_r2, alpha, gain_norm);
    % z_r2_tilde = norm_trans(z_r2, z_r1, z_r2, alpha, gain_norm);
    z_r1_tilde = norm_trans(z_r1, z_r1, z_r2, alpha_gain, minus_alpha_gain, gain_norm);
    z_r2_tilde = norm_trans(z_r2, z_r1, z_r2, alpha_gain, minus_alpha_gain, gain_norm);
    z_r1_tilde = z_r1_tilde(C{k});
    z_r2_tilde = z_r2_tilde(C{k});

    L = {[z_r1_tilde, z_r2_tilde]};  % The first column is the lower bound, and the second column is the upper bound.
    % P = {false(length(C{k}),2)};  % 可以有序储存，但感觉提升不大甚至会提高运行时间，因为新元素进来还是得排序
    % P = {true(length(C{k}),2)};
    P = {[false(length(C{k}),1) true(length(C{k}),1)]};

    log_cut = zeros(p-1+1,100000);
    pointer_uncut = 1;  % 对于每个矩形，指向未执行的切的开头
    pointer_log_cut = 1;  % 指向log_cut的空白列

    log_z_star_tilde = zeros(p-1,100000);
    pointer_log_z_star_tilde = 1;  % 指向log_z_star_tilde的空白列

    % % left_tilde = ((1-alpha)*gain_norm*f(C{k},:)+alpha*gain_norm*mean(f,1))./(z_r2(C{k})-z_r1(C{k}));
    % % right_tilde = z_r1(C{k})*gain_norm./(z_r2(C{k})-z_r1(C{k}));
    % % left_tilde = (1-alpha)*(f(C{k},:)*gain_norm./(z_r2(C{k})-z_r1(C{k}))) + alpha*(mean(f,1)*gain_norm./(z_r2(C{k})-z_r1(C{k})));  % 确保和norm_trans的运算顺序一致
    % % right_tilde = z_r1(C{k})*gain_norm./(z_r2(C{k})-z_r1(C{k}));
    % left_tilde = minus_alpha_gain*(f(C{k},:)*gain_norm./(z_r2(C{k})-z_r1(C{k}))) + alpha_gain*(mean(f,1)*gain_norm./(z_r2(C{k})-z_r1(C{k})));  % 确保和norm_trans的运算顺序一致
    % right_tilde = z_r1(C{k})*gain_norm*(minus_alpha_gain+alpha_gain)./(z_r2(C{k})-z_r1(C{k}));
    H_left_tilde = cell(1, length(C{k}));
    H_mean = mean(cat(3,H{:}),3);
    for i = 1 : length(C{k})
        H_left_tilde{i} = minus_alpha_gain*(H{C{k}(i)}*gain_norm./(z_r2(C{k}(i))-z_r1(C{k}(i)))) + alpha_gain*(H_mean*gain_norm./(z_r2(C{k}(i))-z_r1(C{k}(i))));  % 确保和norm_trans的运算顺序一致
    end
    minus_H_left_tilde = cellfun(@(H) -H, H_left_tilde, 'UniformOutput', false);
    f_left_tilde = minus_alpha_gain*(f(C{k},:)*gain_norm./(z_r2(C{k})-z_r1(C{k}))) + alpha_gain*(mean(f,1)*gain_norm./(z_r2(C{k})-z_r1(C{k})));  % 确保和norm_trans的运算顺序一致
    right_tilde = z_r1(C{k})*gain_norm*(minus_alpha_gain+alpha_gain)./(z_r2(C{k})-z_r1(C{k}));

    % Main loop
    while ~isempty(L) && cputime-start_time_k < runtime0_k-runtime_pre/p

        %{
        % Visualization
        % try close(1); catch; end
        % figure(1)
        clf
        switch p
            case 3
                axis([z_r1_tilde(1) z_r2_tilde(1) z_r1_tilde(2) z_r2_tilde(2)])
                for i = 1 : length(L)
                    patch([L{i}(1,1),L{i}(1,2),L{i}(1,2),L{i}(1,1)],[L{i}(2,1),L{i}(2,1),L{i}(2,2),L{i}(2,2)],rand(1,3),'FaceAlpha',.5,'LineStyle','none')
                end
                try
                    hold on
                    plot(z_cusp_tilde(1),z_cusp_tilde(2),'ko','MarkerSize',10,'LineWidth',1)
                    plot(z_star_tilde(1),z_star_tilde(2),'k*','MarkerSize',10,'LineWidth',1)

                    z_all_tilde = trans_cone(normalize_2(f*Y_P{k}, z_r1, z_r2), alpha);
                    z_all_tilde = z_all_tilde(C{k},:);
                    plot(z_all_tilde(1,:),z_all_tilde(2,:),'r.','MarkerSize',10,'LineWidth',1)
                catch
                end
            case 4
                grid on
                view(135,30)
                axis([z_r1_tilde(1) z_r2_tilde(1) z_r1_tilde(2) z_r2_tilde(2) z_r1_tilde(3) z_r2_tilde(3)])
                xlabel('$f_1$','Interpreter','latex'), ylabel('$f_2$','Interpreter','latex'), zlabel('$f_3$','Interpreter','latex')
                for i = 1 : length(L)
                    pointA = L{i}(:,1);  pointB = L{i}(:,2);
                    vertices = [
                        pointA';
                        [pointA(1), pointB(2), pointA(3)];
                        [pointA(1), pointB(2), pointB(3)];
                        [pointA(1), pointA(2), pointB(3)];
                        pointB';
                        [pointB(1), pointA(2), pointB(3)];
                        [pointB(1), pointA(2), pointA(3)];
                        [pointB(1), pointB(2), pointA(3)];
                        ];
                    faces = [
                        1, 2, 8, 7;  % 底面    pointA逆时针
                        1, 2, 3, 4;  % 左侧面  pointA逆时针
                        5, 6, 4, 3;  % 顶面    pointB逆时针
                        5, 8, 7, 6;  % 右侧面  pointB逆时针
                        5, 3, 2, 8;  % 前面    pointB逆时针
                        1, 7, 6, 4   % 后面    pointA逆时针
                        ];
                    color = rand(1,3);
                    for j = 1:size(faces, 1)
                        patch('Vertices', vertices, 'Faces', faces(j,:), 'FaceColor', color, 'EdgeColor', 'black', 'FaceAlpha',1,'LineStyle','none');
                    end
                end
                try
                    hold on
                    plot3(z_cusp_tilde(1),z_cusp_tilde(2),z_cusp_tilde(3),'ko','MarkerSize',10,'LineWidth',1)
                    plot3(z_star_tilde(1),z_star_tilde(2),z_star_tilde(3),'k*','MarkerSize',10,'LineWidth',1)
                catch
                end
        end
        title(num2str([k length(L) n_model_solved(k) size(Y_P{k},2)]))
        pause(0.5)
        %}
        % [k length(L) n_model_solved(k) size(Y_P{k},2) y_N(k)]
        
        % tmp = sum([P{:}],1);
        % % [sum(tmp(1:2:end)>=p-1)/length(P) sum(tmp(2:2:end)<p-1)/length(P)]
        % sum(tmp(2:2:end)==((1:p-1)'),2)/length(P)  % upper

        % Select the weight
        L_tmp = [];
        while ~isempty(L) && isempty(L_tmp)
            % I_selected = length(L);
            % [~,I_selected] = min(pointer_uncut);
            [~,I_selected] = min(flip(pointer_uncut,1)); I_selected = length(pointer_uncut)-I_selected+1;
            L_tmp = L(I_selected);
            P_tmp = P(I_selected);
            pointer_uncut_tmp = pointer_uncut(I_selected);
            %
            % if ~any(diff(L_tmp{1},[],2)<tol(C{k}))
            % if ~any(diff(L_tmp{1},[],2)<tol_precision)
                % if sum(P_tmp{1}(:,2),1) == 0 || any(diff(L_tmp{1},[],2)<tol_precision)
                if sum(P_tmp{1}(:,2),1) == 0
                    % 删除无效的矩形（待证明）
                    L_tmp = [];
                    P_tmp = [];
                    pointer_uncut_tmp = [];
                else
                    % 检查是否存在与目前的解互不支配的，且目标k的值超过当前y_N(k)的。由于被已知解支配的区域都被删掉了，所以只需检查后者
                    start_time_solver = cputime;
                    % % [~,~,exitflag] = intlinprog(zeros(1,size(f,2)),1:size(f,2), ...  % 只需判断是否可行，和子问题无关
                    % %     [A; -f(k,:); -left_tilde; left_tilde], ...
                    % %     [b; -y_N(k); -L_tmp{1}(:,1)-right_tilde-tol_precision; L_tmp{1}(:,2)+right_tilde-tol_precision], ...  % 不包含等号，这样会更快
                    % %     Aeq,beq,lb,ub,[],options);
                    % % [x0,~,exitflag] = intlinprog(zeros(1,size(f,2)),1:size(f,2), ...  % 只需判断是否可行，和子问题无关
                    % %     [A; -f(k,:); -left_tilde; left_tilde], ...
                    % %     [b; -y_N(k); -L_tmp{1}(:,1)-right_tilde-options_fea.ConstraintTolerance*10; L_tmp{1}(:,2)+right_tilde-options_fea.ConstraintTolerance*10], ...  % 不包含等号，这样会更快    只能用这个，否则会多删了东西。。。
                    % %     Aeq,beq,lb,ub,[],options_fea);
                    % [x0,~,exitflag] = intlinprog(zeros(1,size(f,2)),1:size(f,2), ...  % 只需判断是否可行，和子问题无关
                    %     [A; -f(k,:); -left_tilde; left_tilde], ...
                    %     [b; -y_N(k)-tol; ...  % 不包含等号，且是原始的tol！
                    %     -L_tmp{1}(:,1)-right_tilde; L_tmp{1}(:,2)+right_tilde], ...  % 包含等号
                    %     Aeq,beq,lb,ub,[],options_fea);  % 看起来更加合理
                    Opt = opti('f',zeros(1,size(f,2))', ...
                        'qc', [{-1/2*H{k}}, minus_H_left_tilde, H_left_tilde], [-f(k,:)', -f_left_tilde', f_left_tilde'], [-y_N(k)-tol; -L_tmp{1}(:,1)-right_tilde; L_tmp{1}(:,2)+right_tilde], ...
                        'ineq',A,b, ...
                        'eq',Aeq,beq,'bounds',lb,ub,'xtype',1:size(f,2),'options',options_fea);
                    [x0,~,exitflag] = solve(Opt);

                    runtime_solver = runtime_solver + cputime - start_time_solver;
                    n_model_solved(k) = n_model_solved(k) + 1;
                    switch exitflag
                        case 1
                            % do nothing
                        case 0
                            % 可行检查超时，将其认为具备可行点，此时x0为空
                        case -1
                            L_tmp = [];
                            P_tmp = [];
                            pointer_uncut_tmp = [];
                        case -2
                            error('Unknown error (-2).')
                        otherwise
                            error(['The parameters of the solver are not suitable (' num2str(exitflag) ').'])
                    end
                end
            % end
            %}
            while ~isempty(L_tmp) && ~isempty(pointer_uncut_tmp) && pointer_uncut_tmp(1) <= pointer_log_cut-1
                switch log_cut(end,pointer_uncut_tmp(1))
                    case 1
                        % Rectangular subdivision process (obtained objective vector)
                        [L_tmp, P_tmp] = sub_rects_3(log_cut(1:end-1,pointer_uncut_tmp(1)), z_r2_tilde, L_tmp, P_tmp, 1, true);
                        % [L_tmp, is_remove] = remove_rects_2(L_tmp, log_cut(1:end-1,pointer_uncut_tmp(1)), z_r2_tilde, tol(C{k}));
                        [L_tmp, is_remove] = remove_rects_2(L_tmp, log_cut(1:end-1,pointer_uncut_tmp(1)), z_r2_tilde, 0);  % 只能用这个，否则会多删了东西。。。
                    case 2
                        % Rectangular subdivision process (cusp of contour)
                        [L_tmp, P_tmp] = sub_rects_3(z_r1_tilde, log_cut(1:end-1,pointer_uncut_tmp(1)), L_tmp, P_tmp, 2, false);
                        % [L_tmp, is_remove] = remove_rects_2(L_tmp, z_r1_tilde, log_cut(1:end-1,pointer_uncut_tmp(1)), tol_precision);  % 不能使用tol(C{k})，因为尖点不一定是实际存在的目标向量
                        [L_tmp, is_remove] = remove_rects_2(L_tmp, z_r1_tilde, log_cut(1:end-1,pointer_uncut_tmp(1)), 0);  % 只能用这个，否则会多删东西。。。
                end
                P_tmp(is_remove) = [];
                pointer_uncut_tmp = ones(length(L_tmp),1)*pointer_uncut_tmp(1) + 1;
            end
            L = [L(1:I_selected-1) L_tmp L(I_selected+1:end)];
            P = [P(1:I_selected-1) P_tmp P(I_selected+1:end)];
            pointer_uncut = [pointer_uncut(1:I_selected-1); pointer_uncut_tmp; pointer_uncut(I_selected+1:end)];
        end
        if isempty(L)
            break
        end
        % z_cusp_tilde = L_tmp{1}(:,2);
        z_cusp_tilde = L_tmp{end}(:,2);
        % weight = gain_weight ./ z_cusp_tilde;
        % weight = weight / sum(weight);
        % weight(isnan(weight)) = 1;  % 应该不可能出现0/0的情况
        % Solve the boundary subproblem
        % % A_bd = [[zeros(size(A,1),1) A]; [-ones(p-1,1) weight.*left_tilde]];  % 求解器精度设置使得必须扩大倍数
        % % b_bd = [b; weight.*right_tilde];
        % A_bd = [[zeros(size(A,1),1) A]; [-ones(p-1,1) f_left_tilde./z_cusp_tilde]];  % 直接套入表达式，减少误差
        % b_bd = [b; right_tilde./z_cusp_tilde];
        % x0 = [];
        Q_bd = cell(1,length(H_left_tilde));
        for i = 1 : length(z_cusp_tilde)
            Q_bd{i} = 1/2*[zeros(1,size(H_left_tilde{i},2)+1);zeros(size(H_left_tilde{i},1),1),H_left_tilde{i}./z_cusp_tilde(i)];  % 直接套入表达式，减少误差
        end
        l_bd = [-ones(p-1,1) f_left_tilde./z_cusp_tilde];
        r_bd = right_tilde./z_cusp_tilde;
        A_bd = [zeros(size(A,1),1) A];
        b_bd = b;
        % if ~isempty(x0)
        %     x0 = [max((weight.*left_tilde*x0-weight.*right_tilde)); x0];
        % end
        start_time_solver = cputime;
        % [x_star,~,exitflag] = intlinprog(f_bd,2:length(f_bd),A_bd,b_bd,Aeq_bd,beq_bd,lb_bd,ub_bd,x0,options);
        Opt = opti('f',f_bd, ...  % 无需转置
            'qc', Q_bd, l_bd', r_bd, ...
            'ineq',A_bd,b_bd, ...
            'eq',Aeq_bd,beq_bd,'bounds',lb_bd,ub_bd,'xtype',2:length(f_bd),'options',options);
        [x_star,~,exitflag] = solve(Opt);

        runtime_solver = runtime_solver + cputime - start_time_solver;
        n_model_solved(k) = n_model_solved(k) + 1;
        switch exitflag
            case 1  % Optimal
                % Obtained objective vector
                x_star = round(x_star(2:end));  % 忽略第一个变量；舍入，可能导致解变得不可行，但ismember需要
                % % z_star_tilde = trans_cone(normalize_2(f*x_star, z_r1, z_r2)*gain_norm, alpha)/gain_norm;
                % % z_star_tilde = norm_trans(f*x_star, z_r1, z_r2, alpha, gain_norm);
                % z_star_tilde = norm_trans(f*x_star, z_r1, z_r2, alpha_gain, minus_alpha_gain, gain_norm);
                z_star_tilde = norm_trans(cellfun(@(H)1/2*(x_star')*H*x_star,H)'+f*x_star, z_r1, z_r2, alpha_gain, minus_alpha_gain, gain_norm);
                z_star_tilde = z_star_tilde(C{k});
                % Update the archive and the estimated nadir point
                if ~(~isempty(Y_P{k}) && any(ismember(Y_P{k}', x_star', 'row')))
                    Y_P{k} = [Y_P{k} x_star];
                    f_k_x_star = 1/2*(x_star')*H{k}*x_star + f(k,:)*x_star;
                    if f_k_x_star > y_N(k)
                        y_N(k) = f_k_x_star;
                    end
                    log_cut(:,pointer_log_cut) = [z_star_tilde;1];
                    log_z_star_tilde(:,pointer_log_z_star_tilde) = z_star_tilde;
                    pointer_log_z_star_tilde = pointer_log_z_star_tilde + 1;
                    [L, is_remove] = remove_rects_2(L, z_star_tilde, z_r2_tilde, 0);  % 立刻删除
                    P(is_remove) = []; pointer_uncut(is_remove) = [];
                else  % 该矩形顶点构成的子问题的最优解就是之前的解，也就是说这个矩形的顶点就是方向向量的尖点
                    log_cut(:,pointer_log_cut) = [z_cusp_tilde;2];
                    [L, is_remove] = remove_rects_2(L, z_r1_tilde, z_cusp_tilde, 0);  % 立刻删除
                    P(is_remove) = []; pointer_uncut(is_remove) = [];
                end
                pointer_log_cut = pointer_log_cut + 1;
            case -1
                % ***-2,-9可能可以说明这个矩形有问题
                [L, is_remove] = remove_rects_2(L, L_tmp{end}(:,1), L_tmp{end}(:,2), 0);
                P(is_remove) = []; pointer_uncut(is_remove) = [];
                % if exitflag == -2
                %     error('All solutions are infeasible.')
                % elseif exitflag == -9
                %     error(['exitflag = ' num2str(exitflag) '.'])
                % end
            case -2
                error('The best function value is infinite.')
            otherwise
                error(['The parameters of the solver are not suitable (' num2str(exitflag) ').'])
        end
    end
    runtime_k(k,1) = cputime-start_time_k;
    if (runtime0_k-runtime_pre/p)-(cputime-start_time_k) <= 0
        runtime_k(k,2) = 1;
    end
end
runtime = cputime-start_time;

%%
end
