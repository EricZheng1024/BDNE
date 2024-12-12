classdef ECRNSGAII < ALGORITHM
% <multi/many> <real/binary/permutation>
% NSGA-II with emphasized critical regions strategy

%------------------------------- Reference --------------------------------
% "Nadir point estimation for many-objective optimization problems based on
% emphasized critical regions"
% 
% reproduced by: Ruihao Zheng
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Generate random population
            Population = Problem.Initialization();
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring  = OperatorGA_2(Population(MatingPool),{1,20,1,20});
                [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Problem.N);
            end
        end
    end
end

%%
function [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Population,N)
% The environmental selection

    % Non-dominated sorting
    objs = Population.objs;
    [FrontNo,MaxFNo] = NDSort(objs,Population.cons,N);
    Next = find(FrontNo < MaxFNo);
    
    % Individual number assignment
    N_last_sel = N-length(Next);
    n1 = repmat(floor(N_last_sel / size(objs,2)), size(objs,2), 1);
    tmp = mod(N_last_sel, size(objs,2));
    n1(randperm(size(objs,2),tmp)) = n1(1) + 1;
    
    % Calculate the crowding distance of each solution
    CrowdDis = -inf(1,length(Population));
    for i = 1 : MaxFNo
        index = find(FrontNo == i);
        [~,I] = sort(Population(index).objs, 1);
        index2 = I(end,:);
        if i == MaxFNo
            index_last_sort = I;
            epsilon = sqrt( sum((max(objs(index,:),1)-min(objs(index,:),1)).^2) )' ./ n1;
        end
        [~,I] = sort(I, 1);
        CrowdDis(index) = min(I, [], 2)+1;
        CrowdDis(index(index2)) = 1;
    end
    CrowdDis = -CrowdDis;
    
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);
    
    % epsilon-clearing
    Last = find(FrontNo==MaxFNo);
    for i = 1 : size(objs,2)
        if n1(i) > 0
            counter = 1;
            S = false(1,length(Last));
            S(index_last_sort(end,i)) = true;
            for k = 1 : N_last_sel-1
                if counter >= n1(i)
                    break;
                end
                % if S(index_last_sort(k,i)) == false && max(abs(objs(Last(S),i)-objs(Last(index_last_sort(k,i)),i))) > epsilon(i)
                if S(index_last_sort(k,i)) == false && max(max(abs(objs(Last(S),:)-objs(Last(index_last_sort(k,i)),:)),[],2)) > epsilon(i)
                    S(index_last_sort(k,i)) = true;
                    counter = counter + 1;
                end
            end
            Next = [Next Last(S)];
        end
    end
    
    % Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end
