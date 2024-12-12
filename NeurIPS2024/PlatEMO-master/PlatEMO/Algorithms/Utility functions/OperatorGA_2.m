function Offspring = OperatorGA_2(Parent,Parameter)
%OperatorGA - Crossover and mutation operators of genetic algorithm.
%
%   If the offspring generating by crossover is identical to parents,
%   mutation must be used.

%------------------------------- Reference --------------------------------
% [1] K. Deb, K. Sindhya, and T. Okabe, Self-adaptive simulated binary
% crossover for real-parameter optimization, Proceedings of the Annual
% Conference on Genetic and Evolutionary Computation, 2007, 1187-1194.
% [2] K. Deb and M. Goyal, A combined genetic adaptive search (GeneAS) for
% engineering design, Computer Science and informatics, 1996, 26: 30-45.
% [3] L. Davis, Applying adaptive algorithms to epistatic domains,
% Proceedings of the International Joint Conference on Artificial
% Intelligence, 1985, 162-164.
% [4] D. B. Fogel, An evolutionary approach to the traveling salesman
% problem, Biological Cybernetics, 1988, 60(2): 139-144.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Parameter setting
    if nargin > 1
        [proC,disC,proM,disM] = deal(Parameter{:});
    else
        [proC,disC,proM,disM] = deal(1,20,1,20);
    end
    if isa(Parent(1),'SOLUTION')
        calObj = true;
        Parent = Parent.decs;
    else
        calObj = false;
    end
    Parent1 = Parent(1:floor(end/2),:);
    Parent2 = Parent(floor(end/2)+1:floor(end/2)*2,:);
    [N,D]   = size(Parent1);
    Problem = PROBLEM.Current();
    
    switch Problem.encoding
        case 'real'
            %% Genetic operators for real encoding
            % Simulated binary crossover
            beta = zeros(N,D);
            mu   = rand(N,D);
            beta(mu<=0.5) = (2*mu(mu<=0.5)).^(1/(disC+1));
            beta(mu>0.5)  = (2-2*mu(mu>0.5)).^(-1/(disC+1));
            beta = beta.*(-1).^randi([0,1],N,D);
            beta(rand(N,D)<0.5) = 1;
            beta(repmat(rand(N,1)>proC,1,D)) = 1;
            Offspring = [(Parent1+Parent2)/2+beta.*(Parent1-Parent2)/2
                         (Parent1+Parent2)/2-beta.*(Parent1-Parent2)/2];
            % Polynomial mutation
            Lower = repmat(Problem.lower,2*N,1);
            Upper = repmat(Problem.upper,2*N,1);
            Site  = rand(2*N,D) < proM/D;

            % enhance the mutation of offspring identical to parents    使用该策略仍不能完全避免重复解，因为不同父母对可能产生相同的解
            iden_index = all(Offspring == [Parent1; Parent1], 2) | all(Offspring == [Parent2; Parent2], 2);
            for i = 1 : length(iden_index)
                if iden_index(i)
                    Site(i, randi([1 D])) = true;
                end
            end


            % Offspring = Offspring + randn(size(Offspring)).*(Upper-Lower)*0.05;
            % Offspring = min(max(Offspring,Lower),Upper);


            mu    = rand(2*N,D);
            temp  = Site & mu<=0.5;
            Offspring       = min(max(Offspring,Lower),Upper);
            Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                              (1-(Offspring(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
            temp = Site & mu>0.5; 
            Offspring(temp) = Offspring(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                              (1-(Upper(temp)-Offspring(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
        otherwise
            error('Unsupported representation.')
    end
    if calObj
        Offspring = SOLUTION(Offspring);
    end
end