function [y_I, n_model_solved, runtime_solver, T, y_PT, I_PT] = payoff_table_MOIQP(H,f,A,b,lb,ub,Aeq,beq,p, ...
    options, n_model_solved, runtime_solver)

    % Obtain the ideal objective vector
    x_PT = zeros(size(f,2),p);
    y_I = zeros(p,1);
    for m = 1 : p
        H_m = H{m};
        f_m = f(m,:);
        start_time_solver = cputime;
        % [tmp,~,exitflag] = intlinprog(f_m,1:length(f_m),A,b,Aeq,beq,lb,ub,[],options);
        Opt = opti('qp',H_m,f_m','ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub,'xtype',1:length(f_m),'options',options);
        [tmp,~,exitflag] = solve(Opt);
        runtime_solver = runtime_solver + cputime - start_time_solver;
        n_model_solved(m) = n_model_solved(m) + 1;
        switch exitflag
            case 1
                x_PT(:,m) = round(tmp);
                y_I(m) = 1/2*(x_PT(:,m)')*H{m}*x_PT(:,m) + f(m,:)*x_PT(:,m);  % 采用舍入后的解计算
            case -1
                error(['The best solution of ' num2str(m) '-th objective does not exist.'])
            case -2
                error(['The best value of ' num2str(m) '-th objective is infinite.'])
            otherwise
                error(['The parameters of the solver are not suitable for optimizing the ' num2str(m) '-th objective. (' num2str(exitflag) ').'])
        end
    end
    % x_PT = round(x_PT);  % 舍入，可能导致解变得不可行
    % y_I = sum(f.*x_PT',2);  % 采用舍入后的解计算

    % Obtain the objective vector contributing to the ideal objective vector
    if nargout > 3
        H_sum = sum(cat(3, H{:}), 3);
        f_sum = sum(f,1);
        for m = 1 : p
            H_m = H{m};
            f_m = f(m,:);
            start_time_solver = cputime;
            % [tmp,~,exitflag] = intlinprog(f_sum,1:length(f_m),A,b,[Aeq; f_m],[beq; y_I(m)],lb,ub,[],options);
            Opt = opti('qp',H_sum,f_sum','qc',1/2*H_m,f_m',y_I(m),'ineq',A,b,'eq',Aeq,beq,'bounds',lb,ub,'xtype',1:length(f_m),'options',options);
            [tmp,~,exitflag] = solve(Opt);
            runtime_solver = runtime_solver + cputime - start_time_solver;
            n_model_solved(m) = n_model_solved(m) + 1;
            switch exitflag
                case 1
                    x_PT(:,m) = tmp;
                case -1
                    % if ~all([A*x_PT(:,m)-b<=options.ConstraintTolerance; abs([Aeq; f_m]*x_PT(:,m)-[beq; y_I(m)])<=options.ConstraintTolerance])  % 可行域中只有一个点时，intlinprog似乎失效
                    %     error('Unknown error (-2).')
                    % end
                    if isempty(A)
                        flag_A = true;
                    else
                        flag_A = A*x_PT(:,m)-b<=0;
                    end
                    if isempty(Aeq)
                        flag_Aeq = true;
                    else
                        flag_Aeq = abs(Aeq*x_PT(:,m)-beq)<=0;
                    end
                    if ~all([1/2*(x_PT(:,m)')*H_m*x_PT(:,m)+f(m,:)*x_PT(:,m)-y_I(m)<=0;    flag_A;    flag_Aeq])  % 可行域中只有一个点时，scip也似乎失效
                        error('Unknown error (-2).')
                    end
                case -2
                    error(['The best value of ' num2str(m) '-th objective is infinite.'])
                otherwise
                    error(['The parameters of the solver are not suitable for optimizing the ' num2str(m) '-th objective. (' num2str(exitflag) ').'])
            end
        end
        x_PT = round(x_PT);  % 舍入，可能导致解变得不可行
        T = f*x_PT;  % T(k,m) 第m个解在第k个目标上的值
        for k = 1 : p
            for m = 1 : p
                T(k,m) = T(k,m) + 1/2*(x_PT(:,m)')*H{k}*x_PT(:,m);
            end
        end
        [y_PT, I_PT] = max(T,[],2);
    end
end