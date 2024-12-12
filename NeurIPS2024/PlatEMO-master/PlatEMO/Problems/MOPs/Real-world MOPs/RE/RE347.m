classdef RE347 < PROBLEM
% <multi/many> <real>
% RE3-4-7 (i.e., RE37)

%------------------------------- Reference --------------------------------
% "An Easy-to-use Real-world Multi-objective Optimization Problem Suite"
%--------------------------------------------------------------------------

    properties(Access = private)
        name;
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            obj.M = 3;
            obj.D = 4;
            obj.lower = [0 0 0 0];
            obj.upper = [1 1 1 1];
            obj.encoding = 'real';
            obj.name = 'RE37';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            PopObj = zeros(size(PopDec, 1), obj.M);
            for i = 1 : size(PopDec, 1)
                PopObj(i,:) = feval(obj.name, PopDec(i,:));
            end
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            % R = load(['reference_points_' obj.name '.dat']);  % uniform reference points

            R = load(['reference_points_nadir_RE347.mat']);  % nonuniform reference points but better nadir objective vector estimation
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