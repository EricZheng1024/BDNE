classdef DNPE < ALGORITHM
% <multi/many> <real/binary/permutation>
% Decomposition-based nadir point estimation method
% lambda --- 100 --- 

%------------------------------- Reference --------------------------------
% Y. Sun, G. G. Yen, and Z. Yi, IGD indicator-based evolutionary algorithm
% for many-objective optimization problems, IEEE Transactions on
% Evolutionary Computation, 2019, 23(2): 173-187.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            lambda = Algorithm.ParameterSet(100);

            %% Nadir point estimation
            Population = Problem.Initialization();
            while Algorithm.NotTerminated(Population)
                Offspring  = OperatorGA_2(Population(randi(end,1,Problem.N)),{1,20,1,20});
                Population = [Population,Offspring];
                fit = Fitness(Population.objs,lambda);
                [~,rank]   = sort(fit,1);
                Population = Population(unique(rank(1:ceil(Problem.N/Problem.M),:)));
            end
        end
    end
end

function fit = Fitness(PopObj,lambda)
% Calculate the objective value of each solution on each single-objective
% optimization problem in nadir point estimation

    fit   = zeros(size(PopObj));
    for i = 1 : size(PopObj,2)
        fit(:,i) = abs(PopObj(:,i)) + lambda*sum(PopObj(:,[1:i-1,i+1:end]).^2,2);
    end
end
