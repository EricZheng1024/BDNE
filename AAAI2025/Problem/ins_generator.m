clear

%% general config
% ns = [20,10,10,30,15,15, 200,50,50,240,60,60, 80,40,40,100,50,50, 50,25,25,60,30,30];
% ps = [3,4,5,3,4,5,3,4,5,3,4,5,3,4,5,3,4,5,3,4,5,3,4,5];
% problem_names = {'AP','AP','AP','AP','AP','AP', ...
%     'KP','KP','KP','KP','KP','KP', ...
%     'ILP','ILP','ILP','ILP','ILP','ILP', ...
%     'SPP','SPP','SPP','SPP','SPP','SPP'};

ns = [4,4,4,6,6,6];
ps = [3,4,5,3,4,5];
problem_names = {'IQP','IQP','IQP','IQP','IQP','IQP'};

n_ins = 30;


%%
pathstr = fileparts(mfilename('fullpath'));  % 本m文件所在的路径
cd(pathstr);  % 更改当前活动目录路径
for ii = 1 : length(ns)
n = ns(ii);
p = ps(ii);
problem_name = problem_names{ii};

%%
[~,~] = mkdir(fullfile('mine',problem_name));
for i_ins = 1 : n_ins
    switch problem_name
        case 'KP'
            % f = -randi([10,100],[p,n]);
            % A = randi([10,100],[1,n]);
            f = -randi([1,1000],[p,n]);
            A = randi([1,1000],[1,n]);
            b = sum(A,2)/2;
            lb = zeros(n,1);
            ub = ones(n,1);
            Aeq = [];
            beq = [];
        case 'AP'
            f = randi([1,20],[p,n*n]);
            Aeq = zeros(2*n,n*n);
            for i = 1 : n
                Aeq(i,n*(i-1)+1:n*i) = 1;  % row
                Aeq(n+i,i:n:end) = 1;  % colume
            end
            beq = ones(2*n,1);
            lb = zeros(n*n,1);
            ub = ones(n*n,1);
            A = [];
            b = [];
        case 'ILP'
            % "An exact method for computing the nadir values in multiple objective linear programming" 
            m = round(n / 2);
            % obj
            f = randi([-100,-1],[p,n]);
            prob = rand([p,n]);
            tmp = randi([0,100],[p,n]);
            f(prob>0.2) = tmp(prob>0.2);
            % con
            A = randi([-100,-1],[m,n]);
            prob = rand([m,n]);
            A(prob>0.1&prob<0.2) = 0;
            tmp = randi([1,100],[m,n]);
            A(prob>0.2) = tmp(prob>0.2);
            b = zeros(m,1);
            for i = 1 : m
                b(i) = randi([100,max(sum(A(i,:)),100)],1);
            end
            % rest
            lb = zeros(n,1);
            ub = inf(n,1);
            Aeq = [];
            beq = [];
        case 'SPP'
            % "Approaches for multiobjective combinatorial optimization problems" (2007)
            % obj
            stage = [];
            while n-2-sum(stage) >= round((n-2)*0.12)
                stage = [stage randi(round((n-2)*[0.08,0.12]),1)];
            end
            stage(end+1) = n-2-sum(stage);  % (n_1-1),(n_2-n_1),...,(n_s-n_{s-1})
            stage = cumsum(stage) + 1; stage(end) = stage(end) + 1;  % n_1,...,n_s    稍微作了修改，将第一个点加入了第一个stage，最后一个点加入了最后的stage
            c = cell(1,p);
            for k = 1 : p
                c{k} = (n*100)*ones(n);
                % c{k} = zeros(n);
                for i = 1 : n
                    s_i = find(stage>=i,1);
                    for j = 1 : n
                        s_j = find(stage>=j,1);
                        if s_i == s_j && i < j
                            c{k}(i,j) = randi([10,50],1);
                        elseif s_i+1 == s_j
                            c{k}(i,j) = randi([30,100],1);
                        end
                    end
                end
                c{k} = c{k}(:);
                % c{k}(c{k}==0) = sum(c{k});
            end
            f = [c{:}]';
            % con
            % Aeq = zeros(n*(n-1),n*n);
            % beq = zeros(n*(n-1),1);  % default 0
            % for k = 2 : n
            %     for i = 1 : n
            %         index = (k-1-1)*n+i;
            %         Aeq(index,i:n:end) = 1;  % ij
            %         Aeq(index,(1:n)+(i-1)*n) = Aeq(index,(1:n)+(i-1)*n) - 1;  % ji
            %         if i == 1
            %             beq(index) = 1;
            %         elseif i == k
            %             beq(index) = -1;
            %         end
            %     end
            % end
            Aeq = zeros(n,n*n);
            beq = zeros(n,1);  % default 0
            for i = 1 : n
                Aeq(i,i:n:end) = 1;  % ij
                Aeq(i,(1:n)+(i-1)*n) = Aeq(i,(1:n)+(i-1)*n) - 1;  % ji
                if i == 1
                    beq(i) = 1;
                elseif i == n  % k = n
                    beq(i) = -1;
                end
            end
            % rest
            lb = zeros(n*n,1);
            ub = ones(n*n,1);
            A = [];
            b = [];
        case 'IQP'
            % Modified from the ILP.
            m = round(n / 2);
            % obj  modified
            H = cell(1,p);
            for i = 1 : p
                % H{i} = randi([-50,-1],[n,n])*2;  % "*2"确保目标函数值为整数
                % prob = rand([n,n]);
                % tmp = randi([0,50],[n,n])*2;
                % H{i}(prob>0.2) = tmp(prob>0.2);
                % H{i} = randi([0,50],[n,n])*2;  % "*2"确保目标函数值为整数

                H{i} = randi([0,3],[n,n]);  % 构造半正定矩阵
                H{i} = (H{i}')*H{i}*2;
            end
            % f = randi([-100,-1],[p,n]);
            % prob = rand([p,n]);
            % tmp = randi([0,100],[p,n]);
            % f(prob>0.2) = tmp(prob>0.2);
            f = randi([-20,-1],[p,n]);
            prob = rand([p,n]);
            tmp = randi([0,20],[p,n]);
            f(prob>0.2) = tmp(prob>0.2);
            % con
            A = randi([-100,-1],[m,n]);
            prob = rand([m,n]);
            A(prob>0.1&prob<0.2) = 0;
            tmp = randi([1,100],[m,n]);
            A(prob>0.2) = tmp(prob>0.2);
            b = zeros(m,1);
            for i = 1 : m
                b(i) = randi([100,max(sum(A(i,:)),100)],1);
            end
            % rest  modified
            lb = -5*ones(n,1);
            ub = 5*ones(n,1);
            Aeq = [];
            beq = [];
        otherwise
            error('Undefined problem.')
    end
    switch problem_name
        case {'ILP', 'IQP'}
            filename = fullfile('mine',problem_name,[problem_name '_p-' num2str(p) '_n-' num2str(n) '_m-' num2str(m) '_ins-' num2str(i_ins)]);
        otherwise
            filename = fullfile('mine',problem_name,[problem_name '_p-' num2str(p) '_n-' num2str(n) '_ins-' num2str(i_ins)]);
    end
    if exist([filename '.mat'], 'file')
        disp(['已存在' filename])
        continue
    end
    switch problem_name
        case 'ILP'
            save(fullfile('mine',problem_name,[problem_name '_p-' num2str(p) '_n-' num2str(n) '_m-' num2str(m) '_ins-' num2str(i_ins)]), 'f', 'A', 'b', 'Aeq', 'beq', 'lb', 'ub', 'p', 'n');
        case 'IQP'
            save(fullfile('mine',problem_name,[problem_name '_p-' num2str(p) '_n-' num2str(n) '_m-' num2str(m) '_ins-' num2str(i_ins)]), ...
                'H', 'f', 'A', 'b', 'Aeq', 'beq', 'lb', 'ub', 'p', 'n');
        otherwise
            save(fullfile('mine',problem_name,[problem_name '_p-' num2str(p) '_n-' num2str(n) '_ins-' num2str(i_ins)]), 'f', 'A', 'b', 'Aeq', 'beq', 'lb', 'ub', 'p', 'n');
    end
end

%%
end
