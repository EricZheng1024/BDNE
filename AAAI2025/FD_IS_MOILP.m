% clear, close all

function [y_N,Y_P,n_model_solved,runtime,runtime_solver,runtime_k] = FD_IS_MOILP(f,A,b,lb,ub,Aeq,beq,min_finterval,p, ...
    filename, ...
    runtime0_k,options)

%% Config
% general_config;
tol = min_finterval / 2;

n_model_solved = zeros(p+1, 1);
runtime_solver = 0;
runtime_k = zeros(p+1,2);
start_time = cputime;


%% Ideal objective vector
% [y_I, n_model_solved, runtime_solver] = payoff_table(f,A,b,lb,ub,Aeq,beq,p, ...
%     options, n_model_solved, runtime_solver);
[y_I, n_model_solved, runtime_solver, ~, y_PT, I_PT] = payoff_table(f,A,b,lb,ub,Aeq,beq,p, ...
    options, n_model_solved, runtime_solver);



%% Upper objective values
[y_U, n_model_solved, runtime_solver] = upper_bound(f,A,b,lb,ub,Aeq,beq,p, ...
    options, n_model_solved, runtime_solver);


%% Facet detection
% Init
runtime_pre = cputime-start_time;
M = 100*(y_U-y_I);  % should let the ideal objective vector be zero
Y_M = (M.*eye(p))';  % col vectors
L = Y_M(:);
m_prime = Y_M(:,1);
y_N_tilde = y_I;
Y_E_prime = [];
V = [];
F = [];
Y_P = [];
while ~isempty(L) && cputime-start_time < runtime0_k*p
    % [1 size(L,2) n_model_solved(p+1) size(Y_E_prime,2)]

    % Step (1)    可以有不同的选择方式
    % FIFO
    index_R = 1;
    % LIFO
    % index_R = size(L,2);
    % The third queue discipline in "Approaches for multiobjective combinatorial optimization problems"
    % index_R = find(utility_R(L, Y_E_prime) == nchoosek(p,2)-1, 1);
    % if isempty(index_R)
    %     index_R = 1;  % FIFO
    % end
    % Specify R
    R = reshape(L(:,index_R),p,[]);
    if isempty(V)
        V = R(:);
    else
        V = union(V',R(:)','rows','stable')';
    end
    % Step (2)
    lambda = sym('lambda',[p,1]);
    res = solve([lambda'*R(:,2:end)==lambda'*R(:,1)  sum(lambda)==1]);
    % res = solve([sum(lambda.*R(:,2:end),1)==sum(lambda.*R(:,1),1)  sum(lambda)==1]);
    lambda = zeros(p,1);
    for q = 1 : p
        tmp = eval(res.(['lambda' num2str(q)]));
        if ~isempty(tmp)
            lambda(q) = tmp;
        else
            lambda(q) = 0;
        end
    end

    if all(lambda > 0)
        % Step (3)    MOIP(lambda)
        start_time_solver = cputime;
        [x_star,~,exitflag] = intlinprog(lambda'*f,1:size(f,2),A,b,Aeq,beq,lb,ub,[],options);
        runtime_solver = runtime_solver + cputime - start_time_solver;
        n_model_solved(p+1) = n_model_solved(p+1) + 1;
        switch exitflag
            case {1,3}
                x_star = round(x_star);  % 舍入，可能导致解变得不可行
                r_star_raw = f*x_star;
                r_star = r_star_raw - y_I;  % 平移理想点
                if any(ismember(R', r_star', 'rows'))
                    % Step (4)
                    if isempty(F)
                        F = R(:);
                    else
                        F = union(F',R(:)','rows','stable')';
                    end
                    is_step_6_7 = false;
                else
                    % Step (5)
                    if isempty(Y_E_prime)
                        Y_E_prime = r_star;
                        Y_P = x_star;
                    else
                        Y_E_prime = union(Y_E_prime',r_star','rows','stable')';
                        Y_P = union(Y_P',x_star','rows','stable')';
                    end
                    y_N_tilde = max([y_N_tilde,r_star_raw],[],2);
                    is_step_6_7 = true;
                end
            case -2
                error('All solutions are infeasible.')
            case -3
                error('The best function value is infinite.')
            otherwise
                error(['The parameters of the solver are not suitable (' num2str(exitflag) ').'])
        end
    else
        % if size(Y_E_prime,2) > p
        %     % Step (8)
        %     tmp = setdiff(Y_E_prime',R','rows');
        %     r_star = tmp(1,:)';
        % else
        %     % Step (9)
        %     r_star = Y_M(:,1);
        % end
        if size(Y_E_prime,2) > p
            % Step (8)
            tmp = setdiff(Y_E_prime',R','rows');
            r_star = tmp(1,:)';
        else
            % Step (9)
            % if ~any(ismember(Y_M', m_prime(:,index_R)', 'rows'))  % for debug
            %     error('m_prime is not in Y_M.')
            % end
            r_star = m_prime(:,index_R);
        end
        is_step_6_7 = true;
    end
    % Step (6) & (7)
    if is_step_6_7
        for q = 1 : p
            R_q = R;
            m_prime_tmp = R_q(:,q);  % 不支持LIFO
            R_q(:,q) = r_star;
            if any(ismember(Y_M',R_q','rows')) && ~any(ismember(V',(R_q(:))','rows'))
                tmp = size(L,2);
                L = union(L',R_q(:)','rows','stable')';
                if size(Y_E_prime,2) <= p && tmp ~= size(L,2)
                    if ~any(ismember(Y_M', m_prime_tmp', 'rows'))
                        m_prime_tmp = intersect(R_q',Y_M','rows')'; m_prime_tmp = m_prime_tmp(:,randi(end));
                    end
                    m_prime = [m_prime m_prime_tmp];
                end
            end
        end
    end
    % Step (10)  Step (6) & (7)只是在L后append，index_R指向元素不变
    L(:,index_R) = [];
    if size(Y_E_prime,2) <= p
        m_prime(:,index_R) = [];
    end
end
runtime_k(p+1,1) = cputime-start_time-runtime_pre;
if cputime-start_time >= runtime0_k*p
    runtime_k(p+1,2) = 1;
end


%% Interior search
runtime_pre = cputime-start_time;
y_N = y_N_tilde;
F0 = F;
while ~isempty(F) && cputime-start_time < runtime0_k*p
    % [2 size(F,2) n_model_solved(1:p)']

    % Step (1)
    % index_R = 1;  % FIFO  可以有不同的选择方式
    index_R = size(F,2);  % LIFO
    R = reshape(F(:,index_R),p,[]);
    % Step (2)
    F(:,index_R) = [];
    % Step (3)    MOIP(R)
    % [y_N_tmp,Y_P_tmp,n_model_solved_tmp,~,runtime_solver_tmp] = KL_MOILP(f,[A; f],[b; max(R,[],2)+y_I],lb,ub,Aeq,beq,min_finterval,p, ...
    %     [], ...
    %     runtime0_k,options, ...
    %     [], [], [], y_U);  % 虽然新增约束会改变y_I,I_PT,y_U，但理论上不影响；而y_PT则需要更新，所以连带y_I,I_PT需要更新 
    % [y_N_tmp,Y_P_tmp,n_model_solved_tmp,~,runtime_solver_tmp,runtime_k_tmp] = KL_MOILP(f,[A; f],[b; max(R,[],2)+y_I],lb,ub,Aeq,beq,min_finterval,p, ...
    %     [], ...
    %     runtime0_k*p-(cputime-start_time),options, ...  % 不使用runtime0_k-(cputime-start_time)/p，因为可能存在超时目标减去其他目标的冗余时间后不超总时间
    %     y_I, y_PT, I_PT, y_U);  % 虽然新增约束会改变y_I,y_PT,I_PT,y_U，但理论上不影响 
    [y_N_tmp,Y_P_tmp,n_model_solved_tmp,~,runtime_solver_tmp,runtime_k_tmp] = KL_MOILP(f,[A; f],[b; max(R,[],2)+y_I],lb,ub,Aeq,beq,min_finterval,p, ...
        [], ...
        runtime0_k*p-(cputime-start_time),options, ...  % 不使用runtime0_k-(cputime-start_time)/p，因为可能存在超时目标减去其他目标的冗余时间后不超总时间
        y_I, y_N, I_PT, y_U);  % 虽然新增约束会改变y_I,y_PT,I_PT,y_U，但理论上不影响    用目前的y_N代替y_PT
    y_N = max([y_N,y_N_tmp],[],2);
    tmp = [Y_P_tmp{:}];
    if ~isempty(tmp)
        Y_P = union(Y_P',[Y_P_tmp{:}]','rows','stable')';
    end
    n_model_solved(1:p) = n_model_solved(1:p) + n_model_solved_tmp;
    runtime_solver = runtime_solver + runtime_solver_tmp;
    runtime_k(1:p,1) = runtime_k(1:p,1) + runtime_k_tmp(:,1);
end
if cputime-start_time >= runtime0_k*p
    runtime_k(runtime_k(1:p,1)<=runtime0_k-runtime_pre/p,2) = 1;
end
runtime = cputime-start_time;

%%
end


%%
function res = utility_R(L, Y_E_prime)

    res = zeros(length(L),1);
    if isempty(Y_E_prime)
        return
    end
    for i = 1 : length(L)
        R = L{i};
        for t = 1 : size(R,2)
            r_t = R(:,t);
            for s = 1 : size(R,2)
                r_s = R(:,s);
                tmp = CF_prime(r_s, L, Y_E_prime);
                if ~isempty(tmp)
                    res(i) = res(i) + size(intersect(r_t',tmp','rows','stable'),1);
                end
            end
        end
    end
    res = round(res / 2);
end


function y_prime = CF_prime(y, L, Y_E_prime)
    y_prime = [];
    for i = 1 : length(L)
        if any(ismember(L{i}',y','rows'))
            tmp = intersect(Y_E_prime',L{i}','rows','stable');
            if ~isempty(tmp)
                if isempty(y_prime)
                    y_prime = tmp;
                else
                    y_prime = union(y_prime,tmp,'rows','stable');
                end
            end
        end
    end
    y_prime = y_prime';
end
