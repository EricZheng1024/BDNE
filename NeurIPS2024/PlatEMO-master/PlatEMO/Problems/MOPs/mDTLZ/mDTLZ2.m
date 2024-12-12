classdef mDTLZ2 < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
% Modified DTLZ2

%------------------------------- Reference --------------------------------
% "On Scalable Multiobjective Test Problems With Hardly Dominated
% Boundaries"
%--------------------------------------------------------------------------

    methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            if isempty(obj.D); obj.D = obj.M+7; end
            obj.lower    = zeros(1,obj.D);
            obj.upper    = ones(1,obj.D);
            obj.encoding = 'real';
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            N = size(PopDec,1);
            M = obj.M;
            g = zeros(N,M);
            for i = 1 : N
                for j = 1 : M
                    phi_k = PopDec(i, M+j-1:M : obj.D);
                    g(i,j) = sum((phi_k - 0.5).^2);
                end
            end
            PopObj = (1+g).*(1-fliplr(cumprod([ones(size(g,1),1),cos(PopDec(:,1:M-1)*pi/2)],2)).*[ones(size(g,1),1),sin(PopDec(:,M-1:-1:1)*pi/2)]);
        end
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
            x_grid = UniformPoint(N,obj.M-1,'grid')';
            R = zeros(N,obj.M);
            for i = 1 : size(x_grid,2)
                x = x_grid(:,i);
                X = zeros(obj.M,1);
                X(1) = 1-prod(cos(0.5*pi*x(1:end)));
                X(obj.M) = 1-sin(0.5*pi*x(1));
                for j = 2 : obj.M-1
                    X(j) = 1-prod(cos(0.5*pi*x(1:obj.M-j)))*sin(0.5*pi*x(obj.M-j+1));
                end
                R(i,:) = X';
            end
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            switch obj.M
                case 2
                    R = obj.GetOptimum(100);
                case 3
                    N_sqrt = 20;
                    tmp = obj.GetOptimum(N_sqrt^2);
                    R = cell(1,3);
                    for i = 0 : (N_sqrt^2-1)
                        R{1}(mod(i,N_sqrt)+1,floor(i/N_sqrt)+1) = tmp(i+1,1);
                        R{2}(mod(i,N_sqrt)+1,floor(i/N_sqrt)+1) = tmp(i+1,2);
                        R{3}(mod(i,N_sqrt)+1,floor(i/N_sqrt)+1) = tmp(i+1,3);
                    end
                otherwise
                    R = [];
            end
        end
    end
end