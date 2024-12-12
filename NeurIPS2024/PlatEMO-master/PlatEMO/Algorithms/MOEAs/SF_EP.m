classdef SF_EP < ALGORITHM
% <multi/many> <real/binary/permutation> <constrained/none>
% All-in-One: Heuristic Methods for Estimating the Nadir Objective Vector
% (except ECR-NSGA-II and DNPE)
% type --- 1.2 --- 
% 
% Author: Ruihao Zheng
% Last modified: 12/12/2024

    methods
        function main(Algorithm,Problem)
            %% Initialization
            % Parameter setting
            type = Algorithm.ParameterSet("EP2");
            if isnumeric(type)
                switch type
                    case 1.2
                        type = "SF2";
                    case 1.3
                        type = "SF3";
                    case 1.4
                        type = "SF4";
                    case 2.2
                        type = "EP2";
                    case 2.3
                        type = "EP3_INF";
                    case 2.4
                        type = "EP4";
                    case 2.5
                        type = "EP5";
                    case 3  % EP1
                        type = "LEX_PAYOFF_1";
                end
            end

            % Generate random population
            Population = Problem.Initialization();
            [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Problem.N,type);

            %% Optimization
            while Algorithm.NotTerminated(Population)
                if type == "LEX_PAYOFF_1" && Problem.FE > Problem.maxFE / 2
                    type = "LEX_PAYOFF_2";
                end
                MatingPool = TournamentSelection(2,Problem.N,FrontNo,-CrowdDis);
                Offspring  = OperatorGA_2(Population(MatingPool));
                [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Problem.N,type);
            end
        end
    end
end


%%
function [Population,FrontNo,CrowdDis] = EnvironmentalSelection(Population,N,type)
% The environmental selection of NSGA-II

%------------------------------- Copyright --------------------------------
% Copyright (c) 2022 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    % Non-dominated sorting
    [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
    Next = FrontNo < MaxFNo;
    
    % Calculate the crowding distance of each solution
    CrowdDis = CrowdingDistance(type,Population.objs,FrontNo);
    
    % Select the solutions in the last front based on their crowding distances
    Last     = find(FrontNo==MaxFNo);
    [~,Rank] = sort(CrowdDis(Last),'descend');
    Next(Last(Rank(1:N-sum(Next)))) = true;
    
    % Population for next generation
    Population = Population(Next);
    FrontNo    = FrontNo(Next);
    CrowdDis   = CrowdDis(Next);
end


%%
function CrowdDis = CrowdingDistance(type,PopObj,FrontNo)
% Rewrite

    [N,M] = size(PopObj);
    if nargin < 2
        FrontNo = ones(1,N);
    end
    CrowdDis = zeros(1,N);

    Fronts   = setdiff(unique(FrontNo),inf);
    for f = 1 : length(Fronts)
        Front = find(FrontNo==Fronts(f));
        Fmax  = max(PopObj(Front,:),[],1);
        Fmin  = min(PopObj(Front,:),[],1);
        PopObj_tmp = (PopObj(Front,:)-Fmin)./(Fmax-Fmin);

        switch type
            case 'SF2'
                res = corner_selection(PopObj_tmp);
                CrowdDis(Front(res>=0)) = 2;
                CrowdDis(Front(res<0)) = 1;
            case 'SF3'
                CrowdDis(Front) = I_r_p_all(PopObj_tmp,inf);
            case 'SF4'
                CrowdDis(Front) = I_SDE_plus_all(PopObj_tmp);
            case 'EP2'
                fit = inf(1,size(PopObj_tmp,1));
                j = 1;  % select one first, and then choose based on distance.
                [~,I] = min(PopObj_tmp(:,j));
                index = find(PopObj_tmp(:,j)==PopObj_tmp(I,j));
                index_rest = setdiff(1:size(PopObj_tmp,1),index);
                count = 0;
                fit(index) = count;
                ext_points = PopObj_tmp(index,:);
                for i = index_rest
                    count = count + 1;
                    dist = mean(pdist2(PopObj_tmp,ext_points),2);
                    [~,I] = max(dist);
                    fit(I) = count;
                    ext_points = [ext_points; PopObj_tmp(I,:)];
                end
                CrowdDis(Front) = -fit;
            case 'EP3_INF'
                fit = zeros(size(PopObj_tmp));
                axis_vector = eye(M);
                axis_vector(axis_vector == 0) = 1e-6;
                for i = 1 : M
                    fit(:,i) = vecnorm(PopObj_tmp./axis_vector(i,:), inf, 2);
                end
                [~,I] = sort(fit,1);
                [~,I] = sort(I,1);
                I = size(PopObj_tmp,1)-I+1;
                CrowdDis(Front) = max(I,[],2);
            case 'EP4'
                fit = zeros(size(PopObj_tmp));
                axis_vector = eye(M);
                for i = 1 : M
                    normW = sqrt(sum(axis_vector(i,:).^2, 2));
                    normObjs = sqrt(sum(PopObj_tmp.^2, 2));
                    Cosine = sum(PopObj_tmp.*axis_vector(i,:), 2) ./ normW ./ normObjs;
                    fit(:,i) = normObjs.*sqrt(1-Cosine.^2);
                end
                [~,I] = sort(fit,1);
                [~,I] = sort(I,1);
                I = size(PopObj_tmp,1)-I+1;
                CrowdDis(Front) = max(I,[],2);
            case 'EP5'
                fit = zeros(size(PopObj_tmp));
                axis_vector = eye(M);
                for i = 1 : M
                    fit(:,i) = pdist2(PopObj_tmp,axis_vector(i,:),'cosine');
                end
                [~,I] = sort(fit,1);
                [~,I] = sort(I,1);
                I = size(PopObj_tmp,1)-I+1;
                CrowdDis(Front) = max(I,[],2);
            case {'LEX_PAYOFF_1','LEX_PAYOFF_2'}
                fit = zeros(size(PopObj_tmp));
                switch type
                    case 'LEX_PAYOFF_1'
                        fit = PopObj_tmp;
                    case 'LEX_PAYOFF_2'
                        axis_vector = eye(M)*1e6+1;
                        for i = 1 :M
                            fit(:,i) = sum(PopObj_tmp.*axis_vector(i,:), 2);
                        end
                end
                [~,I] = sort(fit,1);
                [~,I] = sort(I,1);
                I = size(PopObj_tmp,1)-I+1;
                CrowdDis(Front) = max(I,[],2);

            otherwise
                error('Undefined method.')
        end

        % for i = 1 : M
        %     [~,Rank] = sortrows(PopObj(Front,i));
        % 
        % 
        %     CrowdDis(Front(Rank(1)))   = inf;
        %     CrowdDis(Front(Rank(end))) = inf;
        %     for j = 2 : length(Front)-1
        %         PopObj(Front(Rank(j+1)),i)
        %         CrowdDis(Front(Rank(j))) = CrowdDis(Front(Rank(j)))+(PopObj(Front(Rank(j+1)),i)-PopObj(Front(Rank(j-1)),i))/(Fmax(i)-Fmin(i));
        %     end
        % end
    end

end


%% SF2
function res = corner_selection(P)
    res = -1*ones(1, size(P,1));
    z_min = min(P,[],1);
    z_max = max(P,[],1);
    res(any(P<z_min+z_max/5, 2)) = 1;
end


%% SF3
function res = I_r_p_all(P,p)
    res = zeros(1, size(P,1));
    for i = 1 : size(P,1)
        res(i) = I_r_p(P(i,:), setdiff(P,P(i,:),'rows'), p);
    end
end

function res = I_r_p(y, P, p)
    res = zeros(1, size(P,1));
    for i = 1 : size(P,1)
        if any(y<P(i,:))  % y has at least one better objective function value
        % if 1
            res(i) = vecnorm(R_u_v(y,P(i,:)),p,2);
        else
            res(i) = -vecnorm(R_u_v(P(i,:),y),p,2);
        end
    end
    res = min(res);
end

function res = R_u_v(u, v)
    res = max(v./u-1, 0);
end


%% SF4
function res = I_SDE_plus_all(P)
    res = zeros(1,size(P,1));
    for i = 1 : size(P,1)
        res(i) = I_SDE_plus(P(i,:), setdiff(P,P(i,:),'rows'));
    end
end

function res = I_SDE_plus(y, P)
    s_y = WS(y);
    s_P = WS(P);
    res = nan(1,size(P,1));
    % res = inf(1,size(P,1));
    for i = 1 : size(P,1)
        if s_y - s_P(i) > 1e-9  % WS_p > WS_q
            res(i) = pdist2(y, obj_shift(P(i,:),y));
        end

        % res(i) = pdist2(y, obj_shift(P(i,:),y));  % I_SDE
    end

    % if all(isnan(res))
    %     for i = 1 : size(P,1)
    %         res(i) = pdist2(y, obj_shift(P(i,:),y));  % I_SDE
    %     end
    % end
    
    res = min(res);
end

function res = WS(objs)
    res = sum(objs,2);
end

function res = obj_shift(q, p)
    res = max(q, p);
end
