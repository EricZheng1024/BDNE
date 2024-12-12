classdef CRE531 < PROBLEM
% <multi/many> <real>
% CRE5-3-1 (i.e., CRE51), the modified version of RE6-3-1

%------------------------------- Reference --------------------------------
% "An Easy-to-use Real-world Multi-objective Optimization Problem Suite"
%--------------------------------------------------------------------------

    properties(Access = private)
        name;
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 5;
            obj.D = 3;
            obj.lower = [0.01 0.01 0.01];
            obj.upper = [0.45 0.10 0.10];
            obj.encoding = 'real';
            obj.name = 'CRE51';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = zeros(size(PopDec, 1), obj.M);
            for i = 1 : size(PopDec, 1)
                [PopObj(i,:),PopCon] = feval(obj.name, PopDec(i,:));
                PopObj(i,:) = PopObj(i,:) + 1e6*sum(PopCon);  % penalty
            end
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            % R = load(['reference_points_RE61.dat']);  % uniform reference points
            % R = R(R(:,obj.M+1)<=0,1:obj.M);

            R = load(['reference_points_nadir_CRE531.mat']);  % nonuniform reference points but better nadir objective vector estimation
            R = R.R;
        end
        %% Generate the image of Pareto front
        function DrawObj(obj,Population)
            % rewrite
            ax = Draw(Population.objs,{'\it f\rm_1','\it f\rm_2','\it f\rm_3'});
            if obj.M == 2
                plot(ax,obj.optimum(:,1),obj.optimum(:,2),'.k');
            elseif obj.M == 3
                plot3(ax,obj.optimum(:,1),obj.optimum(:,2),obj.optimum(:,3),'.k');
            end
        end
    end
end