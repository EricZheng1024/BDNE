function [y_U, n_model_solved, runtime_solver] = upper_bound(f,A,b,lb,ub,Aeq,beq,p, ...
    options, n_model_solved, runtime_solver)

    y_U = zeros(p,1);
    for m = 1 : p
        f_m = f(m,:);
        start_time_solver = cputime;
        [tmp,~,exitflag] = intlinprog(-f_m,1:length(f_m),A,b,Aeq,beq,lb,ub,[],options);
        runtime_solver = runtime_solver + cputime - start_time_solver;
        n_model_solved(m) = n_model_solved(m) + 1;
        switch exitflag
            case 1
                tmp = round(tmp);  % 舍入，可能导致解变得不可行
                y_U(m) = f(m,:)*tmp;  % 采用舍入后的解计算
            case -2
                error(['The worst solution of ' num2str(m) '-th objective does not exist.'])
            case -3
                error(['The worst value of ' num2str(m) '-th objective is infinite.'])
            otherwise
                error(['The parameters of the solver are not suitable for maximizing the ' num2str(m) '-th objective. (' num2str(exitflag) ').'])
        end
    end
end