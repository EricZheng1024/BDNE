function [y_I, n_model_solved, runtime_solver, T, y_PT, I_PT] = payoff_table(f,A,b,lb,ub,Aeq,beq,p, ...
    options, n_model_solved, runtime_solver)

    % Obtain the ideal objective vector
    x_PT = zeros(size(f,2),p);
    for m = 1 : p
        f_m = f(m,:);
        start_time_solver = cputime;
        [tmp,~,exitflag] = intlinprog(f_m,1:length(f_m),A,b,Aeq,beq,lb,ub,[],options);
        runtime_solver = runtime_solver + cputime - start_time_solver;
        n_model_solved(m) = n_model_solved(m) + 1;
        switch exitflag
            case 1
                x_PT(:,m) = tmp;
            case -2
                error(['The best solution of ' num2str(m) '-th objective does not exist.'])
            case -3
                error(['The best value of ' num2str(m) '-th objective is infinite.'])
            otherwise
                error(['The parameters of the solver are not suitable for optimizing the ' num2str(m) '-th objective. (' num2str(exitflag) ').'])
        end
    end
    x_PT = round(x_PT);  % 舍入，可能导致解变得不可行
    y_I = sum(f.*x_PT',2);  % 采用舍入后的解计算

    % Obtain the objective vector contributing to the ideal objective vector
    if nargout > 3
        f_sum = sum(f,1);
        for m = 1 : p
            f_m = f(m,:);
            start_time_solver = cputime;
            [tmp,~,exitflag] = intlinprog(f_sum,1:length(f_m),A,b,[Aeq; f_m],[beq; y_I(m)],lb,ub,[],options);
            runtime_solver = runtime_solver + cputime - start_time_solver;
            n_model_solved(m) = n_model_solved(m) + 1;
            switch exitflag
                case 1
                    x_PT(:,m) = tmp;
                case -2
                    if ~all([A*x_PT(:,m)-b<=options.ConstraintTolerance; abs([Aeq; f_m]*x_PT(:,m)-[beq; y_I(m)])<=options.ConstraintTolerance])  % 可行域中只有一个点时，intlinprog似乎失效
                        error('Unknown error (-2).')
                    end
                case -3
                    error(['The best value of ' num2str(m) '-th objective is infinite.'])
                otherwise
                    error(['The parameters of the solver are not suitable for optimizing the ' num2str(m) '-th objective. (' num2str(exitflag) ').'])
            end
        end
        x_PT = round(x_PT);  % 舍入，可能导致解变得不可行
        T = f*x_PT;  % T(k,m) 第m个解在第k个目标上的值
        [y_PT, I_PT] = max(T,[],2);
    end
end