% clear, close all

function [y_N,Y_P,n_model_solved,runtime,runtime_solver,runtime_k] = NPDA_MOILP(f,A,b,lb,ub,Aeq,beq,min_finterval,p, ...
    filename, ...
    runtime0_k,options)

%% Config
% general_config;
tol = min_finterval/2;  % P,Q问题的不等式约束应该是不带等号的，因此减去一个极小值
n_model_solved = zeros(p, 1);
runtime_solver = 0;
runtime_k = zeros(p,2);
start_time = cputime;


%% Pay-off table
[y_I, n_model_solved, runtime_solver, T, y_PT, I_PT] = payoff_table(f,A,b,lb,ub,Aeq,beq,p, ...
    options, n_model_solved, runtime_solver);


%% Upper objective values
[y_U, n_model_solved, runtime_solver] = upper_bound(f,A,b,lb,ub,Aeq,beq,p, ...
    options, n_model_solved, runtime_solver);
y_U = y_U + 1;  % Arbitrary upper bounds, delta=1. It does not need to be exact.


%% Main procedure
runtime_pre = cputime-start_time;
y_N = y_PT;
Y_P = cell(1, p);
f_sum = sum(f,1);
for k = 1 : p
    % Init
    start_time_k = cputime;
    m = I_PT(k);
    % m = mod(k,p)+1;
    C = setdiff(1:p,[m,k]);  % 从小到大
    y_I_bar = y_I(C);
    y_U_bar = y_U(C);
    phi = T(C,m);  % 关于C的索引，下同
    C_prime = find(y_I_bar<phi);
    L = cell(1, (2^length(C_prime))-1);  % 数量是杨辉三角行的和、不包含第一列
    count = 1;
    for i = 1 : length(C_prime)  % No empty set is selected since i>=1.
        Js = nchoosek(C_prime, i);
        for j = 1 : size(Js,1)
            J = Js(j,:);
            J_d = setdiff(1:length(C),J);
            R_J = zeros(p-2,2);  % The first column is the lower bound, and the second column is the upper bound.
            R_J(J,:) = [y_I_bar(J) phi(J)];
            R_J(J_d,:) = [phi(J_d) y_U_bar(J_d)];
            L{count} = R_J;
            count = count + 1;
        end
    end
    V = cal_rect_vols(L);  % 可以有序储存，但感觉提升不大甚至会提高运行时间，因为新元素进来还是得排序

    % Main loop
    f_m = f(m,:);
    while ~isempty(L) && cputime-start_time_k < runtime0_k-runtime_pre/p
        % Visualization
        %{
        try close(1); catch; end
        figure(1)
        switch p
            case 3
                for i = 1 : length(L)
                    patch([L{i}(1),L{i}(1),L{i}(2),L{i}(2)],[0 1 1 0],rand(1,3),'FaceAlpha',.5,'LineStyle','none')
                end
                axis([y_I(C(1)) y_U(C(1)) y_I(C(2)) y_U(C(2))])
                try
                    hold on
                    plot(f_bar_x_star(1),f_bar_x_star(2),'ko','MarkerSize',10,'LineWidth',1)
                    plot(u_i(1),u_i(2),'k*','MarkerSize',10,'LineWidth',1)
                catch
                end
            case 4
                for i = 1 : length(L)
                    patch([L{i}(1,1),L{i}(1,2),L{i}(1,2),L{i}(1,1)],[L{i}(2,1),L{i}(2,1),L{i}(2,2),L{i}(2,2)],rand(1,3),'FaceAlpha',.5,'LineStyle','none')
                end
                axis([y_I(C(1)) y_U(C(1)) y_I(C(2)) y_U(C(2))])
                try
                    hold on
                    plot(f_bar_x_star(1),f_bar_x_star(2),'ko','MarkerSize',10,'LineWidth',1)
                    plot(u_i(1),u_i(2),'k*','MarkerSize',10,'LineWidth',1)
                catch
                end
            case 5
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
                grid on
                view(135,30)
                axis([y_I(C(1)) y_U(C(1)) y_I(C(2)) y_U(C(2)) y_I(C(3)) y_U(C(3))])
                xlabel('$f_1$','Interpreter','latex'), ylabel('$f_2$','Interpreter','latex'), zlabel('$f_3$','Interpreter','latex')
                try
                    hold on
                    plot3(f_bar_x_star(1),f_bar_x_star(2),f_bar_x_star(3),'ko','MarkerSize',10,'LineWidth',1)
                    plot3(u_i(1),u_i(2),u_i(3),'k*','MarkerSize',10,'LineWidth',1)
                catch
                end
        end
        title(['Number of rectangles (' num2str(k) '-th objective): ' num2str(length(L))])
        pause(0.5)
        %}
        % [k length(L) n_model_solved(k) size(Y_P{k},2)]


        [~,I_highV] = max(V);
        u_i = L{I_highV}(:,2);

        % P
        A_P = [A; f(k,:); f(C,:)];
        b_P = [b; y_U(k); u_i-tol];
        start_time_solver = cputime;
        [tmp,~,exitflag] = intlinprog(f_m,1:length(f_m),A_P,b_P,Aeq,beq,lb,ub,[],options);
        runtime_solver = runtime_solver + cputime - start_time_solver;
        n_model_solved(k) = n_model_solved(k) + 1;
        switch exitflag
            case {1,3}  % Optimal
                tmp = round(tmp);
                z_star = f_m*tmp;  % 如果不准确，会使得R没有可行解
                % R
                A_R = [A; f(k,:); f(C,:)];
                b_R = [b; y_U(k); u_i-tol];
                Aeq_R = [Aeq; f_m];
                beq_R = [beq; z_star];
                start_time_solver = cputime;
                [x_star,~,exitflag] = intlinprog(f_sum,1:length(f_m),A_R,b_R,Aeq_R,beq_R,lb,ub,[],options);
                runtime_solver = runtime_solver + cputime - start_time_solver;
                n_model_solved(k) = n_model_solved(k) + 1;
                switch exitflag
                    case {1,3}
                        x_star = round(x_star);  % 舍入，可能导致解变得不可行，但ismember需要
                    case -2
                        if all([A_R*tmp-b_R<=options.ConstraintTolerance; abs(Aeq_R*tmp-beq_R)<=options.ConstraintTolerance])  % 可行域中只有一个点时，intlinprog似乎失效
                            x_star = tmp;
                        else
                            error('Unknown error (-2).')
                        end
                    case -3
                        error('Unknown error (-3).')
                    otherwise
                        error(['The parameters of the solver are not suitable for solving R (' num2str(exitflag) ').'])
                end
                f_bar_x_star = f(C,:)*x_star;
                if ~(~isempty(Y_P{k}) && any(ismember(Y_P{k}', x_star', 'row')))
                    % Update the estimated nadir point
                    Y_P{k} = [Y_P{k} x_star];
                    f_k_x_star = f(k,:)*x_star;
                    if f_k_x_star > y_N(k)
                        y_N(k) = f_k_x_star;
                    end
                    % Rectangular subdivision process
                    % 
                end
                [L, V] = sub_rects_2(f_bar_x_star, u_i, L, V);  % 因为u_i会变，所以放这里
                [L, V] = remove_rects(L, V, f_bar_x_star, u_i, 1e-6);
            case -2  % All solutions are infeasible
                [L, V] = remove_rects(L, V, y_I_bar, u_i, 1e-6);
            case -3
                error('The best function value is infinite.')
            otherwise
                error(['The parameters of the solver are not suitable for solving P (' num2str(exitflag) ').'])
        end
    end
    runtime_k(k,1) = cputime-start_time_k;
    if (runtime0_k-runtime_pre/p)-(cputime-start_time_k) <= 0
        runtime_k(k,2) = 1;
    end
end
runtime = cputime-start_time;

% cd(fileparts(mfilename('fullpath')));  % 更改当前活动目录路径
% [~,~] = mkdir('results_NPDA');
% [~,tmp,~] = fileparts(filename);
% save(fullfile('results_NPDA',['NPDA_' tmp]), "y_N", "Y_P", "n_model_solved", "runtime");

%%
end
