% clear, close all

function [y_N,Y_P,n_model_solved,runtime,runtime_solver,runtime_k] = KL_MOILP(f,A,b,lb,ub,Aeq,beq,min_finterval,p, ...
    filename, ...
    runtime0_k,options, ...
    y_I, y_PT, I_PT, y_U)

%% Config
% general_config;
tol = min_finterval/2;

n_model_solved = zeros(p, 1);
runtime_solver = 0;
runtime_k = zeros(p,2);
start_time = cputime;


%% Pay-off table & Upper objective values
if ~exist('y_I','var') || isempty(y_I)
    [y_I, n_model_solved, runtime_solver, ~, y_PT, I_PT] = payoff_table(f,A,b,lb,ub,Aeq,beq,p, ...
        options, n_model_solved, runtime_solver);  % z_n^{IP(n)} is y_I; z_n^{IP(m)} is y_PT(m)
end

if ~exist('y_U','var') || isempty(y_U)
    [y_U, n_model_solved, runtime_solver] = upper_bound(f,A,b,lb,ub,Aeq,beq,p, ...
        options, n_model_solved, runtime_solver);
end


%% Main procedure
% 与原文记号保持一致且以最大化问题形式求解，因为子问题函数还未确定转为最小化时应具备的形式
runtime_pre = cputime-start_time;
g_star = 0;
M = 0.5/options.IntegerTolerance;  % M*options.IntegerTolerance要小于1，否则会使得求解出错
gain = 1e4;  % 由于epsilon较小，所以需要增益（与M的设置不冲突）
f = -f;  % 将MOP转换为最大化求解器求解的形式
y_U_0 = y_U;  % make positive integer
y_I = -y_I + y_U_0;
y_PT = -y_PT + y_U_0;
y_U = -y_U + y_U_0;
% y_N = y_I;  % 不能使用y_PT，因为构成其的目标向量未必Pareto最优
y_N = y_PT;
Y_P = cell(1, p);
for n = 1 : p
    % Init
    start_time_k = cputime;
    m = I_PT(n);
    if m == n
        m = mod(n,p)+1;
    end
    J = setdiff(1:p,[m,n]);  % I
    t = 0;
    % lz = y_PT(J) - 1;  % lz_j^k (p-2 rows, index j; t cols, index k) -> lz    for minimization
    lz = y_PT(J) + 1;  % lz_j^k (p-2 rows, index j; t cols, index k) -> lz
    uz_n = y_PT(n);  % uz_n
    lz_n = y_U(n);  % lz_n
    epsilon = gain / sum(y_I([J n]));
    epsilon_n_P = gain / (y_I(n)*sum(y_I([J n])));
    epsilon_n_D = gain / (y_I(n)*(y_I(m)+y_I(n)));
    epsilon_m = gain / (y_I(m)+y_I(n));
    flag = 0;  % 0 1 2

    % Main loop
    while (uz_n-lz_n) / (y_I(n)-lz_n) > g_star  &&  flag ~= 2  &&  cputime-start_time_k < runtime0_k-runtime_pre/p
        if flag ~= -1
            % uz_n_star = floor((uz_n+g_star*y_I(n))/(1-g_star)) + 1;  % uz_n^{star}    for minimization
            uz_n_star = ceil((uz_n-g_star*y_I(n))/(1-g_star)) - 1;  % uz_n^{star}
            flag = 0;
            while flag == 0
                t = t + 1;
                n_slack = t*(p-2);
                % Constraints
                A_P = [A zeros(size(A,1),n_slack);     repmat(-f(J,:),[t 1]) M*eye(n_slack);   f(n,:) zeros(1,n_slack); -f(n,:) zeros(1,n_slack)];
                b_P = [b;     M-lz(:)+repmat(y_U_0(J),[t 1]);      uz_n_star-y_U_0(n); -lz_n+y_U_0(n)];
                tmp = zeros(t,n_slack);  % \sum_{j\neq m,n} y_{jk}    y_{jk} -> y(:)
                for i = 1 : t
                    tmp(i,(1:(p-2))+(i-1)*(p-2)) = 1;
                end
                Aeq_P = [Aeq zeros(size(Aeq,1),n_slack); zeros(t,size(f,2))  tmp];
                beq_P = [beq; ones(t,1)];
                lb_P = [lb; zeros(n_slack,1)];
                ub_P = [ub; ones(n_slack,1)];
                % P
                options.MaxTime = max((runtime0_k-runtime_pre/p)-(cputime-start_time_k), 1);
                start_time_solver = cputime;
                [x_star,~,exitflag] = intlinprog(  -[gain*f(m,:)+epsilon_n_P*f(n,:)+epsilon*sum(f(J,:),1) zeros(1,n_slack)], ...  % 子问题需要最大化，remark: MOP已转换为最大化求解器求解的形式
                    1:(size(f,2)+n_slack), ...
                    A_P,b_P,Aeq_P,beq_P,lb_P,ub_P,[],options);
                runtime_solver = runtime_solver + cputime - start_time_solver;
                n_model_solved(n) = n_model_solved(n) + 1;
                switch exitflag
                    % case 1  % Optimal
                    case {1,3}  % Optimal  等式约束可能不绝对满足，可能是IntegerTolerance的设置
                        x_star = round(x_star(1:size(f,2)));  % 舍入，可能导致解变得不可行
                        dz_t = f*x_star;
                        A_D = [A; -f];
                        b_D = [b; -dz_t];
                        Aeq_D = Aeq;
                        beq_D = beq;
                        lb_D = lb;
                        ub_D = ub;
                        % D
                        options.MaxTime = max((runtime0_k-runtime_pre/p)-(cputime-start_time_k), 1);  % 一旦之前设置了，这里也需要重新设置！
                        start_time_solver = cputime;
                        [x_check,~,exitflag] = intlinprog(  -(gain*sum(f(J,:),1)+epsilon_m*f(m,:)+epsilon_n_D*f(n,:)), ...
                            1:size(f,2), ...
                            A_D,b_D,Aeq_D,beq_D,lb_D,ub_D,[],options);
                        runtime_solver = runtime_solver + cputime - start_time_solver;
                        n_model_solved(n) = n_model_solved(n) + 1;
                        switch exitflag
                            case {1,3}  % Optimal
                                x_check = round(x_check);  % 舍入，可能导致解变得不可行
                                z_t = f*x_check;
                                lz = [lz z_t(J)+y_U_0(J)+1];
                                % if isequal(dz_t,z_t)  % 判断目标值而不是解本身
                                if all(abs(dz_t-z_t)<tol)
                                    uz_n = z_t(n) + y_U_0(n);
                                    flag = 1;
                                end
                                if ~(~isempty(Y_P{n}) && any(ismember(Y_P{n}', x_check', 'row')))
                                    Y_P{n} = [Y_P{n} x_check];
                                end
                            case {0,2}
                                % if (runtime0_k-runtime_pre/p)-(cputime-start_time_k) <= 0  % 模型求解时间上限要大于该目标允许运行时间上限
                                %     flag = -1;
                                % end
                                flag = -1;  % 不搞条件判断，防止求解器停止条件0,2时，仍然有极少时间上的裕量。配合(if flag~=-1)，可以让其一直浪费时间直到裕量被消耗完。
                                % break
                            case {-2,-3}
                                error(['Unknown error. (' num2str(exitflag) ').'])
                            otherwise
                                error(['The parameters of the solver are not suitable for solving D (' num2str(exitflag) ').'])
                        end
                    case -2  % All solutions are infeasible
                        % lz_n = uz_n_star - 1;  % lz_n    for minimization
                        lz_n = uz_n_star + 1;  % lz_n
                        flag = 2;
                    case {0,2}
                        % if (runtime0_k-runtime_pre/p)-(cputime-start_time_k) <= 0  % 模型求解时间上限要大于该目标允许运行时间上限
                        %     flag = -1;
                        % end
                        flag = -1;
                        % break
                    otherwise
                        error(['The parameters of the solver are not suitable for solving P (' num2str(exitflag) ').'])
                end
                % Update the estimated nadir point
                % y_N(n) = (lz_n + uz_n) / 2;
                y_N(n) = uz_n;

                % [n n_model_solved(n) size(Y_P{n},2) -(y_N(n)-y_U_0(n))]
            end
        end
    end
    runtime_k(n,1) = cputime-start_time_k;
    if (runtime0_k-runtime_pre/p)-(cputime-start_time_k) <= 0
        runtime_k(n,2) = 1;
    end
end
y_N = -(y_N - y_U_0);
y_I = -(y_I - y_U_0);
y_PT = -(y_PT - y_U_0);
y_U = -(y_U - y_U_0);

runtime = cputime-start_time;

%%
end
