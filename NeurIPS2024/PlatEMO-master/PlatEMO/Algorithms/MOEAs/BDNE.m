classdef BDNE < ALGORITHM
% <multi/many> <real/binary/permutation>
% Boundary decomposition for nadir objective vector estimation
% mu --- 100 --- bound of trade-off
% gen_pass_max --- 200 --- iteration number of LLOPs
% lambda_def_ULOP ---  --- 
% T --- 0 --- 
% 
% Author: Ruihao Zheng
% Last modified: 12/12/2024

    methods
        function main(Algorithm,Problem)
            %% Initialization
            % Parameter setting
            [mu, gen_pass_max, lambda_def_ULOP, T] = Algorithm.ParameterSet(100, 200, 4+floor(3*log(Problem.M-2)), 0);
            alpha = Problem.M / ((Problem.M-1) * mu + Problem.M);
            if T == 0  % let Problem.N be the total number of LLOPs
                lambda_def_ULOP = floor(Problem.N/Problem.M)-1;
                Problem.N = (lambda_def_ULOP+1)*Problem.M;
            else
                Problem.N = (lambda_def_ULOP+1)*Problem.M*T;
            end
            N_subp = (lambda_def_ULOP+1)*Problem.M;
            N_nad = lambda_def_ULOP*Problem.M;

            % Initial the CMA models of the ULOP
            weight_lower = zeros(1,Problem.M-2);
            weight_upper = ones(1,Problem.M-2);
            C = cov([eye(Problem.M-2); ones(1,Problem.M-2)/(Problem.M-1); zeros(1,Problem.M-2)]);
            [~,D] = eig(C);
            Sigma_ULOP = struct('s',num2cell(1:Problem.M),'lambda',lambda_def_ULOP,'m',1/(Problem.M-1)*ones(1,Problem.M-2), ...
                'sigma',0.3,'sigma_0',0.3, ...
                'C',C,'diagD_0',diag(D), ...
                'p_c',0,'p_sigma',0,'gen',0,'gen_TolFun',0,'sol_best',[]);
            exitflag_ULOP = zeros(Problem.M,1);

            % Generate random population
            Population = Problem.Initialization();

            % Initial perception
            gen_pass = 1;
            [~,FrontNo,CrowdDis] = EnvironmentalSelection_ECNSGAII(Population,Problem.N);
            while gen_pass < gen_pass_max
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring  = OperatorGA_2(Population(MatingPool));
                [Population,FrontNo,CrowdDis] = EnvironmentalSelection_ECNSGAII([Population,Offspring],Problem.N);
                gen_pass = gen_pass + 1;
            end

            % Determine reference points
            z_r1 = min(Population.objs, [], 1);
            z_r2 = max(Population.objs, [], 1);

            % Generate the boundary weight vectors
            W = zeros(N_nad, Problem.M);
            flag_inj_ULOP = false(N_nad, 1);
            index_wi = cell(1,Problem.M);
            index_wi_2 = zeros(N_nad, 1);
            for i = 1 : Problem.M
                tmp = mvnrnd(Sigma_ULOP(i).m,Sigma_ULOP(i).sigma^2*Sigma_ULOP(i).C,Sigma_ULOP(i).lambda-1);
                [tmp,tmp2] = repair_ULOP(tmp,weight_lower,weight_upper);
                tmp = [tmp 1-sum(tmp,2); ones(1,Problem.M-1)/(Problem.M-1)];
                tmp = [tmp(:,1:i-1) zeros(size(tmp,1),1) tmp(:,i:end)];
                index = (1:lambda_def_ULOP)+(i-1)*lambda_def_ULOP;
                W(index, :) = tmp;
                flag_inj_ULOP(index) = [tmp2; true];
                index_wi{i} = index;
                index_wi_2(index) = i;
            end

            % Generate the auxiliary weight vectors and determine the neighbors
            W_aux = eye(Problem.M);

            % Calculate fitness for reproduction
            [~,Population_bs_opt,Population_aux_opt,sol_subp_ranking] = ...
                EnvironmentalSelection(Population, Problem.N, z_r1, z_r2, alpha, W, W_aux);
            Archive = [Population_bs_opt Population_aux_opt];

            % Signal
            gen_pass = 1;
            flag_final_gen = false;

            % Log
            log_gen_pass = [];


            %% Optimization
            while Algorithm.NotTerminated(Archive)
                % Reproduction
                MatingPool = zeros(1, 2*N_subp);
                for i = 1 : N_subp
                    MatingPool(i:N_subp:end) = TournamentSelection(2,2,sol_subp_ranking(:,i));
                end
                Offspring = OperatorGAhalf_2(Population(MatingPool));

                % Determine reference point
                z_r1 = min([z_r1; Offspring.objs], [], 1);

                % Selection
                [Population,Population_bs_opt,Population_aux_opt,sol_subp_ranking] = ...
                    EnvironmentalSelection([Population,Offspring], Problem.N, z_r1, z_r2, alpha, W, W_aux);
                Archive = [Population_bs_opt Population_aux_opt];

                % Update the weight vector
                gen_pass = gen_pass + 1;
                gen_pass_mean = mean([log_gen_pass, gen_pass]);
                gen_remain = (Problem.maxFE-Problem.FE)/N_subp;
                if ~flag_final_gen && ...
                        ( gen_pass >= gen_pass_max || ...
                        gen_remain < gen_pass_mean )  % *extensible design  the remaining generations of LLOP is smaller than the average generations of LLOP; the first generation of LLOP is too long
                    % Reset
                    log_gen_pass = [log_gen_pass gen_pass];
                    gen_pass = 1;
                    % Calculate subproblem fitness
                    subp_fitness = zeros(N_nad, 1);
                    for i = 1 : Problem.M
                        objs = Population_bs_opt(index_wi{i}).objs;
                        subp_fitness(index_wi{i}) = -objs(:,i);
                    end
                    % Determine the update method
                    if (Problem.maxFE-Problem.FE)/(2*N_subp) >= mean(log_gen_pass)  % not last generation of ULOP
                        % Determine reference point
                        z_r2 = max(Population_bs_opt.objs, [], 1);
                        % Update the CMA models of the ULOP
                        for i = 1 : Problem.M
                            % Penalize the flat landscape
                            W_prime = W(index_wi{i},setdiff(1:Problem.M,i));  % M-1 dimension
                            subp_fitness_i = subp_fitness(index_wi{i});
                            % if length(unique(subp_fitness(index_wi{i}),'stable')) < length(subp_fitness(index_wi{i}))
                                % Map the solutions to boundary weight vectors
                                objs = Population_bs_opt(index_wi{i}).objs;
                                objs_trans = trans_cone(normalize_2(objs, z_r1, z_r2, 0), alpha, 1);
                                W_map = 1./objs_trans(:,setdiff(1:Problem.M,i));
                                W_map = W_map ./ sum(W_map,2);
                                % Deal with exception: objs_trans is equal to 0 (only one non-dominated point)
                                index_nan = any(isnan(W_map),2);
                                if any(index_nan)
                                    W_map(index_nan,:) = repmat(ones(1,Problem.M-1)/Problem.M, sum(index_nan), 1);
                                end
                                % Penalize
                                subp_fitness_i = subp_fitness_i + 1e-6*vecnorm(W_map-W_prime,2,2);  % Penalty term coefficient should be samll enough
                            % end
                            % Update the CMA model
                            [~,ranking] = sort(subp_fitness_i);
                            index_free = 1:Problem.M-2;
                            if exitflag_ULOP(i) ~= -1
                                [Sigma_tmp, exitflag_ULOP(i)] = UpdateCMA(Sigma_ULOP(i),W_prime(ranking,index_free),lambda_def_ULOP,flag_inj_ULOP(index_wi{i}(ranking)));
                            end
                            switch exitflag_ULOP(i)
                                case 0
                                    Sigma_ULOP(i) = Sigma_tmp;
                                case -1
                                    % do nothing
                                case 1
                                    Sigma_ULOP(i) = Sigma_tmp;
                                    Sigma_ULOP(i).s = i;
                                    Sigma_ULOP(i).m = 1/(Problem.M-1)*ones(1,Problem.M-1);
                                    Sigma_ULOP(i).sigma = 0.3;
                                    Sigma_ULOP(i).sigma_0 = Sigma_ULOP(i).sigma;
                                    Sigma_ULOP(i).C = eye(Problem.M-1);
                                    [~,D] = eig(Sigma_ULOP(i).C);
                                    Sigma_ULOP(i).diagD_0 = diag(D);
                            end
                            % Genreate new boundary weight vectors
                            n_new_ULOP = ceil(lambda_def_ULOP/2);
                            n_old_ULOP = lambda_def_ULOP - n_new_ULOP;  % number of retained LLOPs
                            tmp = mvnrnd(Sigma_ULOP(i).m,Sigma_ULOP(i).sigma^2*Sigma_ULOP(i).C,n_new_ULOP);
                            [tmp,tmp2] = repair_ULOP(tmp,weight_lower,weight_upper);
                            tmp = [tmp 1-sum(tmp,2)];
                            W(index_wi{i},setdiff(1:Problem.M,i)) = [W_prime(ranking(1:n_old_ULOP),:); tmp];
                            W(index_wi{i}(1),setdiff(1:Problem.M,i)) = W_map(ranking(1),:);  % map
                            flag_inj_ULOP(index_wi{i}) = [true(n_old_ULOP,1); tmp2];
                        end
                    else
                        for i = 1 : Problem.M
                            [~, I_min_subpf] = min(subp_fitness(index_wi{i}));
                            W(index_wi{i},:) = repmat(W(index_wi{i}(I_min_subpf),:), length(index_wi{i}), 1);
                            % W_aux(i,:) = W(index_wi{i}(I_min_subpf),:);
                        end
                        flag_final_gen = true;
                    end
                end
            end
        end
    end
end


%%
function [Population,FrontNo,CrowdDis] = EnvironmentalSelection_ECNSGAII(Population,N)
% The environmental selection of ECNSGA-II

    % Preprocessing
    Population = Population(randperm(length(Population)));  % shuffle
    objs = Population.objs;
    [~,ia] = unique(objs,"rows","stable");  % eliminate duplicate solutions
    Population = Population(ia);

    % Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
    Next = FrontNo < MaxFNo;
    
    % Calculate the crowding distance of each solution
    CrowdDis = -inf(1,length(Population));
    for i = 1 : MaxFNo
        index = FrontNo == i;
        [~,I] = sort(Population(index).objs, 1);
        [~,I] = sort(I, 1);
        I = max(I,sum(index)-I+1);
        CrowdDis(index) = max(I,[],2);
    end
    
    % Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    % Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end


function [Population,Population_bs_opt,Population_aux_opt,sol_subp_ranking] = ...
    EnvironmentalSelection(Population, N, z_r1, z_r2, alpha, W, W_aux)
% The environmental selection

    % Preprocessing
    Population = Population(randperm(length(Population)));  % shuffle
    objs = Population.objs;
    [objs,ia] = unique(objs,"rows","stable");  % eliminate duplicate solutions
    Population = Population(ia);
    [n, m] = size(objs);
    N_nad = size(W,1);

    % Normalization and transformation
    objs_n = normalize_2(objs, z_r1, z_r2, 0);
    objs_trans = trans_cone(objs_n, alpha, 1);

    % Calculate the subproblem function values
    g = zeros(n, size(W,1));
    for i = 1 : size(W,1)
        g(:,i) = g_tch(objs_trans, zeros(1,m), W(i,:));
    end
    g_aux = zeros(n, m);
    for i = 1 : m
        g_aux(:,i) = g_tch(objs_trans, zeros(1,m), W_aux(i,:));
    end
    g = [g g_aux];

    % Calculate the fitenss
    [~,subp_sol_ranking] = sort(g,1);
    [~,sol_subp_ranking] = sort(subp_sol_ranking,1);
    fitness = min(sol_subp_ranking, [], 2);
    Population_bs_opt = Population(subp_sol_ranking(1, 1:N_nad));
    Population_aux_opt = Population(subp_sol_ranking(1, N_nad+1:end));

    % Population for next generation
    [~,Next] = sort(fitness, 'ascend');
    Next = Next(1 : N);
    Population = Population(Next);
    sol_subp_ranking = sol_subp_ranking(Next,:);
end


%%
function objs_n = normalize_2(objs, zmin, zmax, epsilon)
    objs_n = (objs-zmin) ./ (zmax-zmin) - epsilon;
end


function objs = trans_cone(objs,alpha,p)
    objs = (1-alpha)*objs.^p + alpha*mean(objs.^p,2);
end


function g = g_tch(objs, z, W)
    g = max(abs(objs-z) .* W, [], 2);  % Tchebycheff method
end


%%
function [weight_part,flag_inj] = repair_ULOP(weight_part,weight_lower,weight_upper)
    % Frame constraints
    flag_inj = any(weight_part < weight_lower, 2) | any(weight_part > weight_upper, 2);
    weight_part = min(max(weight_part,weight_lower),weight_upper);
    % sum=1
    tmp = sum(weight_part,2);
    index = tmp > 1;
    if any(index)
        flag_inj(index) = true;
        weight_part(index,:) = weight_part(index,:) ./ tmp(index);
    end
end


%%
function [Sigma, exitflag] = UpdateCMA(Sigma,X,lambda_def,flag_inj)
%Update the CMA model

    n = size(X,2);
    Sigma.gen = Sigma.gen + 1;

    % Calculate the CMA parameters
    % Selection
    mu              = floor(Sigma.lambda/2);
    w               = log((Sigma.lambda+1)/2) - log(1:mu);
    w               = w./sum(w);
    mu_eff          = (sum(w)^2) / sum(w.^2);
    c_m             = 1;
    % Injection
    cy              = sqrt(n) + 2*n/(n+2);  
    delta_max_sigma = 1;
    % Adaptation
    c_c             = (4+mu_eff/n) / (n+4+2*mu_eff/n);
    c_sigma         = (mu_eff+2) / (n+mu_eff+5);
    c_1             = 2 / ((n+1.3)^2+mu_eff);
    c_mu            = min(1-c_1,2*(mu_eff-2+1/mu_eff) / ((n+2)^2+mu_eff));
    d_sigma         = 1 + 2*max(0,sqrt((mu_eff-1)/(n+1))-1) + c_sigma;
    ENI             = sqrt(n)*(1-1/4/n+1/21/n^2);
    
    % Update the CMA model
    % Selection and recombination
    y               = (X(1:mu,:)-Sigma.m)/Sigma.sigma;  % row vector
    flag_inj        = flag_inj(1:mu);
    % y(flag_inj,:)   = min(1, cy./norm(Sigma.C^(-1/2)*y(flag_inj,:)')) * y(flag_inj,:);  % injection
    y(flag_inj,:)   = min(1, cy./vecnorm(Sigma.C^(-1/2)*y(flag_inj,:)',2,1))' .* y(flag_inj,:);  % injection
    y_w             = w*y;
    Sigma.m         = Sigma.m + c_m*Sigma.sigma*y_w;
    % Step-size control
    Sigma.p_sigma   = (1-c_sigma)*Sigma.p_sigma + sqrt(c_sigma*(2-c_sigma)*mu_eff)*Sigma.C^(-1/2)*y_w';
    Sigma.sigma     = Sigma.sigma*exp(min(delta_max_sigma,c_sigma/d_sigma*(norm(Sigma.p_sigma)/ENI-1)));  % injection
    % Covariance matrix adaptation
    h_sigma         = norm(Sigma.p_sigma)./sqrt(1-(1-c_sigma).^(2*Sigma.gen)) < (1.4+2/(n+1))*ENI;
    Sigma.p_c       = (1-c_c)*Sigma.p_c + h_sigma*sqrt(c_c*(2-c_c)*mu_eff)*y_w;
    w_circ          = w;
    Sigma.C         = (1+c_1*((1-h_sigma)*c_c*(2-c_c))-c_1-c_mu*sum(w,2))*Sigma.C + c_1*Sigma.p_c'*Sigma.p_c + c_mu*y'*diag(w_circ)*y;
    Sigma.C         = triu(Sigma.C) + triu(Sigma.C,1)'; % Enforce symmetry

    % Reset the CMA model if possible
    exitflag = 0;
    [B,D] = eig(Sigma.C);
    diagD = diag(D);
    diagC = diag(Sigma.C);
    % Routine
    % NoEffectAxis  = all(Sigma.m==Sigma.m+0.1*Sigma.sigma*sqrt(diagD(mod(Sigma.gen,n)+1))*B(mod(Sigma.gen,n)+1,:));
    NoEffectAxis  = all(Sigma.m==Sigma.m+0.1*Sigma.sigma*sqrt(diagD).*B, 'all');
    % NoEffectCoord = all(Sigma.m==Sigma.m+0.2*Sigma.sigma*diagC');
    NoEffectCoord = all(Sigma.m==Sigma.m+0.2*Sigma.sigma*sqrt(diagC)');
    % TolFun        = Sigma.gen_TolFun >= 10 + ceil(30*n/Sigma.lambda);
    TolX          = all(Sigma.sigma*sqrt(diagC) < Sigma.sigma_0*1e-12) && all(Sigma.sigma*Sigma.p_c < Sigma.sigma_0*1e-12);
    % TolX          = all(Sigma.sigma*sqrt(diagC) < Sigma.sigma_0*1e-6) && all(Sigma.sigma*Sigma.p_c < Sigma.sigma_0*1e-6);
    % Exception
    % ConditionCov    = cond(Sigma.C) > 1e14;
    ConditionCov    = false;
    TolXUp          = any(Sigma.sigma*sqrt(diagD) > 1e4*Sigma.sigma_0*sqrt(Sigma.diagD_0));
    ConditionNan    = isnan(Sigma.lambda) || any(isnan(Sigma.m)) || isnan(Sigma.sigma) || any(isnan(Sigma.C),'all') || any(isnan(Sigma.p_c)) || any(isnan(Sigma.p_sigma));

    % if NoEffectAxis || NoEffectCoord || (TolFun && TolX)
    if NoEffectAxis || NoEffectCoord || TolX
        exitflag = -1;
    elseif ConditionCov || TolXUp || ConditionNan
        exitflag = 1;
    end
    if exitflag ~= 0
        Sigma = struct('s',[],'lambda',lambda_def,'m',[],'sigma',[],'sigma_0',[],'C',[],'diagD_0',[],'p_c',0,'p_sigma',0,'gen',0,'gen_TolFun',0,'sol_best',[]);
    end
end
