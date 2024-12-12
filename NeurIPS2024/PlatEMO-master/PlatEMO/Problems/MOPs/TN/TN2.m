classdef TN2 < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
% MOP based on mDTLZ1 with unsupported non-dominated vectors
% point ---  --- the critical point of the m-th objective. e.g. [0.6;0.8;1.2] for the 3-objective case.

%------------------------------- Reference --------------------------------
% 
%--------------------------------------------------------------------------

    properties(Access = private)
        point;
        points_plane;
        i_pmax;
        coeff_plane;
        coeff_conplane;
        sat_con_side;
    end
    
    methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = (2+1)*obj.M; end
            [obj.point] = obj.ParameterSet([0.9; 0.7; 1.5; 0.9*ones(obj.M-3,1)]);
            obj.lower    = 0*ones(1,obj.D);
            obj.upper    = 1*ones(1,obj.D);
            obj.encoding = 'real';
            
            if obj.M < 3
                error('Undefined objective dimension.');
            end
            
            i_pmax_ = find(obj.point > 1);
            if length(i_pmax_) > 1
                error('Undefined setting of the point.');
            end

            inv_simplex = 1 - eye(obj.M);  % row vector
            i_pmax_left = i_pmax_ - 1;
            if i_pmax_left < 1
                i_pmax_left = obj.M-i_pmax_+1;
            end
            i_pmax_right = i_pmax_ + 1;
            if i_pmax_right > obj.M
                i_pmax_right = i_pmax_ - obj.M + 1;
            end
            rest = inv_simplex(setdiff(1:obj.M,[i_pmax_left,i_pmax_right]),:);
            points_plane_{1} = [rest; obj.point'; inv_simplex(i_pmax_left,:)];
            points_plane_{2} = [rest; obj.point'; inv_simplex(i_pmax_right,:)];

            obj.i_pmax = i_pmax_;
            obj.points_plane = points_plane_;
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            N = size(PopDec,1);
            M = obj.M;
            g = zeros(N,M);
            for i = 1 : N
                for j = 1 : M
                    phi_k = PopDec(i, M+j : M : obj.D);
                    g(i,j) = 10*sum((phi_k - 0.3).^2);
                end
            end

            h = zeros(N,M);
            index_plane1 = find(PopDec(:,obj.M)<=0.5);
            index_plane2 = setdiff(1:N, index_plane1);

            PopDec(index_plane1,obj.M) = PopDec(index_plane1,obj.M)*2;
            PopDec(index_plane2,obj.M) = (1-PopDec(index_plane2,obj.M))*2;

            h(index_plane1,:) = PopDec(index_plane1,1:obj.M) ./ sum(PopDec(index_plane1,1:obj.M),2);
            index_nan = any(isnan(h(index_plane1,:)),2);
            h(index_plane1(index_nan),:) = 1/obj.M;
            h(index_plane1,:) = h(index_plane1,:) * obj.points_plane{1};
            
            h(index_plane2,:) = PopDec(index_plane2,1:obj.M) ./ sum(PopDec(index_plane2,1:obj.M),2);
            index_nan = any(isnan(h(index_plane2,:)),2);
            h(index_plane2(index_nan),:) = 1/obj.M;
            h(index_plane2,:) = h(index_plane2,:) * obj.points_plane{2};
            
            PopObj = (1+g).*h;
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            N_sampling = round(N / 2);
            [theta_plane,N_sampling] = UniformPoint(N_sampling,obj.M);
            theta_plane(theta_plane==1e-6) = 0;
            R = zeros(N_sampling*length(obj.points_plane),obj.M);
            for i = 1 : length(obj.points_plane)
                R(((i-1)*N_sampling+1):i*N_sampling,:) = theta_plane*obj.points_plane{i};
            end
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            switch obj.M
                case 3
                    R = cell(1,obj.M);
                    for i = 1 : obj.M
                        R{i} = [obj.points_plane{1}(:,i)  obj.points_plane{2}(:,i)];  % X, Y, Z
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