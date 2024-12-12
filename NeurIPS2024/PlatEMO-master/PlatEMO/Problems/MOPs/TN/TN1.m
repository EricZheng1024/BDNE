classdef TN1 < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
% MOP constructed from endpoints of a hyperplane.
% points --- 4 --- scale or matrix. For scalar, the set of predefined points. For matrix, a row is a point; The number of points is m or m+1.

%------------------------------- Reference --------------------------------
% 
%--------------------------------------------------------------------------

    properties(Access = private)
        points;
        coeff_plane;
        coeff_conplane;
        sat_con_side;
    end
    
    methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = (2+1)*obj.M; end
            [obj.points] = obj.ParameterSet(4);
            if length(obj.points) == 1
                switch obj.points
                    case 1
                        obj.points = eye(obj.M);
                    case 2
                        obj.points = 1 - eye(obj.M);
                    case 3
                        obj.points = 1 - eye(obj.M);
                        obj.points(1,2) = 2;
                    case 4
                        obj.points = 1 - eye(obj.M);
                        obj.points(1,2) = 1.2;
                        obj.points(2,1) = 0.8;
                        obj.points(2,3) = 0.8;
                    case 5
                        obj.points = [1 0.6 0; 0.1 1 1; 0.9 0.3 0.9];
                        obj.M = 3;
                    otherwise
                        error('Undefined set of points.');
                end
            end
            if obj.M ~= size(obj.points,1)
                error('Objective dimension does not match the number of points.')
            end
            obj.lower    = 0*ones(1,obj.D);
            obj.upper    = 1*ones(1,obj.D);
            obj.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            N = size(PopDec,1);
            M = obj.M;
            g = zeros(N,M);
            for i = 1 : N
                for j = 1 : M
                    phi_k = PopDec(i, M+j: M : obj.D);
                    g(i,j) = 10*sum((phi_k - 0.3).^2);
                end
            end
            h = PopDec(:,1:obj.M) ./ sum(PopDec(:,1:obj.M),2);
            index_nan = any(isnan(h),2);
            h(index_nan,:) = 1/obj.M;
            h = h * obj.points;
            PopObj = (1+g).*h;
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            theta_plane = UniformPoint(N,obj.M);
            theta_plane(theta_plane==1e-6) = 0;
            R = theta_plane*obj.points;
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            switch obj.M
                case 3
                    R = cell(1,obj.M);
                    for i = 1 : obj.M
                        R{i} = [obj.points(:,i)];  % X, Y, Z
                    end
                otherwise
                    R = [];
            end
        end
        function DrawObj(obj,Population)
            % rewrite
            ax = Draw(Population.objs,{'\it f\rm_1','\it f\rm_2','\it f\rm_3'});
            switch obj.M
                case 3
                    % patch(ax, obj.PF{1},obj.PF{2},obj.PF{3},'r','FaceColor', '#9bad6e', 'FaceAlpha', 1, 'EdgeColor', '#6d6e71', 'LineStyle','-','LineWidth',1.2)
                    patch(ax, obj.PF{1},obj.PF{2},obj.PF{3},'k', 'FaceAlpha', 0.5, 'EdgeColor', '#6d6e71', 'LineStyle','-','LineWidth',1.2)
                otherwise
                    % do nothing
            end
        end
    end
end