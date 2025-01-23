clear, close all

pathstr = fileparts(mfilename('fullpath'));  % 本m文件所在的路径
cd(pathstr);  % 更改当前活动目录路径

alg_names = {'NPDA','KL','FD_IS','BDNE'};
ns = [20,10,10,30,15,15, 200,50,50,240,60,60, 80,40,40,100,50,50, 50,25,25,60,30,30];
ps = [3,4,5,3,4,5,3,4,5,3,4,5,3,4,5,3,4,5,3,4,5,3,4,5];
problem_names = {'AP','AP','AP','AP','AP','AP', ...
    'KP','KP','KP','KP','KP','KP', ...
    'ILP','ILP','ILP','ILP','ILP','ILP', ...
    'SPP','SPP','SPP','SPP','SPP','SPP'};

% alg_names = {'NPDA','BDNE'};
% ns = [4,4,4,6,6,6];
% ps = [3,4,5,3,4,5];
% problem_names = {'IQP','IQP','IQP','IQP','IQP','IQP'};

% alg_names = {'BDNE_P2','BDNE_woFC','BDNE'};
% ns = [20,10,10,30,15,15, 200,50,50,240,60,60, 80,40,40,100,50,50, 50,25,25,60,30,30];
% ps = [3,4,5,3,4,5,3,4,5,3,4,5,3,4,5,3,4,5,3,4,5,3,4,5];
% problem_names = {'AP','AP','AP','AP','AP','AP', ...
%     'KP','KP','KP','KP','KP','KP', ...
%     'ILP','ILP','ILP','ILP','ILP','ILP', ...
%     'SPP','SPP','SPP','SPP','SPP','SPP'};


n_ins = 30;
runtime0_k = 3600;


%%
runtimes = zeros(length(problem_names),n_ins,length(alg_names));
runtimes_ = zeros(length(problem_names),n_ins,length(alg_names));
y_N_all = cell(length(problem_names)*n_ins,length(alg_names));
y_N_all_ = cell(length(problem_names)*n_ins,length(alg_names));
filenames = strings(length(problem_names)*n_ins*length(alg_names), 1); count = 0;
for k1 = 1 : length(problem_names)
    for i_ins = 1 : n_ins
        switch problem_names{k1}
            case {'ILP','IQP'}
                ins_name = [problem_names{k1} '_p-' num2str(ps(k1)) '_n-' num2str(ns(k1)) '_m-' num2str(round(ns(k1)/2)) '_ins-' num2str(i_ins)];
            otherwise
                ins_name = [problem_names{k1} '_p-' num2str(ps(k1)) '_n-' num2str(ns(k1)) '_ins-' num2str(i_ins)];
        end
        for k2 = 1 : length(alg_names)
            tmp = runtime0_k*ps(k1);
            try
                filename = fullfile('results', ['results_' alg_names{k2}], [alg_names{k2} '_' ins_name]);
                if k2 == 1
                    count = count + 1;
                    filenames(count) = filename;
                end
                load(filename);
                runtime_ = runtime;
                y_N_ = y_N;
                if runtime > tmp || any(runtime_k(:,2))
                    % runtime = tmp;
                    runtime = inf;
                    runtime_k(logical(runtime_k(:,2)),1)=runtime0_k;  runtime_=sum(runtime_k(:,1));
                    y_N = [];
                end
            catch
                disp(['Missing: ' filename '. Consider it a timeout.'])
                % runtime = tmp;
                runtime = inf;
                y_N = [];
                y_N_ = y_N;
                runtime_ = nan;
            end
            runtimes(k1,i_ins,k2) = runtime;
            runtimes_(k1,i_ins,k2) = runtime_;
            y_N_all{(k1-1)*n_ins+i_ins,k2} = y_N;
            y_N_all_{(k1-1)*n_ins+i_ins,k2} = y_N_;
        end
    end
end


%% 检查未能精确估计的实例
%{
check_y_N = false(size(y_N_all,1),1);
for i = 1 : size(y_N_all,1)
    tmp = false(size(y_N_all{i,end}));
    flag = false;
    if ~isempty(tmp)
        for k2 = 1 : size(y_N_all,2)-1
            if ~isempty(y_N_all{i,k2})
                tmp = tmp | (y_N_all{i,k2}==y_N_all{i,end});  % 只要与其中一个精确算法的结果对上就行
                flag = true;
            end
        end
        if flag  % 存在baseline没超时
            check_y_N(i) = all(tmp);
        else  % baseline都超时
            check_y_N(i) = true;
        end
    else  % BDNC都超时
        check_y_N(i) = true;
    end
end
index_inacc = find(~check_y_N);
y_N_inacc = y_N_all(index_inacc,:);
%}


%% 加载部分结果
% % tmp = runtimes(4:9,1:2,1)';
% tmp = runtimes([4 5 7 8],1:3,1)';
% tmp = round(tmp(:));


%% Statistical result
%
med_runtimes = squeeze(median(runtimes,2));
[~,ranking_med_runtimes] = sort(med_runtimes,2);
[~,ranking_med_runtimes] = sort(ranking_med_runtimes,2);

best_runtimes = squeeze(min(runtimes,[],2));
[~,ranking_best_runtimes] = sort(best_runtimes,2);
[~,ranking_best_runtimes] = sort(ranking_best_runtimes,2);

worst_runtimes = squeeze(max(runtimes,[],2));
[~,ranking_worst_runtimes] = sort(worst_runtimes,2);
[~,ranking_worst_runtimes] = sort(ranking_worst_runtimes,2);

worst_best_runtimes = worst_runtimes-best_runtimes;
[~,ranking_worst_best_runtimes] = sort(worst_best_runtimes,2);
[~,ranking_worst_best_runtimes] = sort(ranking_worst_best_runtimes,2);

std_runtimes = squeeze(std(runtimes_,0,2));
[~,ranking_std_runtimes] = sort(std_runtimes,2);
[~,ranking_std_runtimes] = sort(ranking_std_runtimes,2);

MAD_runtimes = squeeze(mad(runtimes_,1,2));
[~,ranking_MAD_runtimes] = sort(MAD_runtimes,2);
[~,ranking_MAD_runtimes] = sort(ranking_MAD_runtimes,2);

rs_runtimes = zeros(length(problem_names),length(alg_names)-1);
for k1 = 1 : length(problem_names)
    for k2 = 1 : length(alg_names)-1
        [~,rs_runtimes(k1,k2)] = ranksum(runtimes_(k1,:,k2),runtimes_(k1,:,length(alg_names)));
    end
end
%}


%% Table
%
n_digits = 6;
n_dec = 1;


% str_stat = 'best';
% res = string(zeros(length(problem_names),length(alg_names)));
% for k1 = 1 : length(problem_names)
%     for k2 = 1 : length(alg_names)
%         if eval([str_stat,'_runtimes(k1,k2)']) == inf
%             res(k1,k2) = ['\ ', '(', num2str(length(alg_names)), ')'];
%         else
%             % res(k1,k2) = strcat(num2str(eval([str_stat,'_runtimes(k1,k2)']), n_digits), '(', num2str(eval(['ranking_',str_stat,'_runtimes(k1,k2)'])), ')');
%             res(k1,k2) = strcat(num2str(round(eval([str_stat,'_runtimes(k1,k2)']), 0), n_digits), '(', num2str(round(eval(['ranking_',str_stat,'_runtimes(k1,k2)'])), 0), ')');
%             % res(k1,k2) = strcat(num2str(round(eval([str_stat,'_runtimes(k1,k2)']), 1), n_digits), '(', num2str(round(eval(['ranking_',str_stat,'_runtimes(k1,k2)'])), 1), ')');
%         end
%     end
% end


% strs_stat = {'best','worst'};
% res = string(nan(length(problem_names),length(alg_names)*length(strs_stat)));
% for i = 1 : length(strs_stat)
%     str_stat = strs_stat{i};
%     for k1 = 1 : length(problem_names)
%         for k2 = 1 : length(alg_names)
%             if eval([str_stat,'_runtimes(k1,k2)']) == inf
%                 res(k1,length(strs_stat)*(k2-1)+i) = ['TO ', '(', num2str(length(alg_names)), ')'];
%             else
%                 res(k1,length(strs_stat)*(k2-1)+i) = strcat(num2str(round(eval([str_stat,'_runtimes(k1,k2)']), n_dec), n_digits), ...
%                     '(', num2str(round(eval(['ranking_',str_stat,'_runtimes(k1,k2)']), n_dec), n_digits), ')');
%             end
%         end
%     end
% end


% dmed,timeout
% res = string(nan(length(problem_names)+2,length(alg_names)*2));
% for k2 = 1 : length(alg_names)
%     win_eq_lose = zeros(1,3);
%     for k1 = 1 : length(problem_names)
%         if k2 == length(alg_names)
%             tmp = '';
%         else
%             switch rs_runtimes(k1,k2)
%                 case 1
%                     if med_runtimes(k1,k2) > med_runtimes(k1,length(alg_names))
%                         tmp = '-';
%                         win_eq_lose(3) = win_eq_lose(3) + 1;
%                     elseif med_runtimes(k1,k2) < med_runtimes(k1,length(alg_names))
%                         tmp = '+';
%                         win_eq_lose(1) = win_eq_lose(1) + 1;
%                     else
%                         error('what?')
%                     end
%                 case 0
%                     tmp = '=';
%                     win_eq_lose(2) = win_eq_lose(2) + 1;
%             end
%         end
%         if med_runtimes(k1,k2) == inf
%             res(k1,2*(k2-1)+1) = ['\ ', '(', num2str(length(alg_names)), ')', tmp];
%         else
%             if k2 == length(alg_names)
%                 res(k1,2*(k2-1)+1) = strcat(num2str(round(med_runtimes(k1,k2), n_dec), n_digits), ...
%                     '(', num2str(ranking_med_runtimes(k1,k2), 0), ')',tmp);
%             else
%                 res(k1,2*(k2-1)+1) = strcat(num2str(round(med_runtimes(k1,k2)-med_runtimes(k1,length(alg_names)), n_dec), n_digits), ...
%                     '(', num2str(ranking_med_runtimes(k1,k2), 0), ')',tmp);
%             end
%         end
%         res(k1,2*(k2-1)+2) = sum(runtimes(k1,:,k2)==inf);
%     end
%     res(length(problem_names)+1,2*(k2-1)+1) = join(string(win_eq_lose), '/');
%     res(length(problem_names)+1,2*(k2-1)+2) = '';
%     res(length(problem_names)+2,2*(k2-1)+1) = round(mean(ranking_med_runtimes(:,k2)), n_dec);
%     res(length(problem_names)+2,2*(k2-1)+2) = '';
% end


% med,dmed,timeout
res = string(nan(length(problem_names)+2,length(alg_names)*3-1));
for k2 = 1 : length(alg_names)
    win_eq_lose = zeros(1,3);
    for k1 = 1 : length(problem_names)
        if k2 == length(alg_names)
            tmp = '';
        else
            switch rs_runtimes(k1,k2)
                case 1
                    if med_runtimes(k1,k2) > med_runtimes(k1,length(alg_names))
                        tmp = '-';
                        win_eq_lose(3) = win_eq_lose(3) + 1;
                    elseif med_runtimes(k1,k2) < med_runtimes(k1,length(alg_names))
                        tmp = '+';
                        win_eq_lose(1) = win_eq_lose(1) + 1;
                    else
                        error('what?')
                    end
                case 0
                    tmp = '=';
                    win_eq_lose(2) = win_eq_lose(2) + 1;
            end
        end
        if med_runtimes(k1,k2) == inf
            res(k1,3*(k2-1)+1) = ['TO ', '(', num2str(length(alg_names)), ')', tmp];
            res(k1,3*(k2-1)+2) = 'TO ';
        else
            if k2 == length(alg_names)
                res(k1,3*(k2-1)+1) = strcat(num2str(round(med_runtimes(k1,k2), n_dec), n_digits), ...
                    '(', num2str(ranking_med_runtimes(k1,k2), 0), ')',tmp);
            else
                res(k1,3*(k2-1)+1) = strcat(num2str(round(med_runtimes(k1,k2), n_dec), n_digits), ...
                    '(', num2str(ranking_med_runtimes(k1,k2), 0), ')',tmp);
                res(k1,3*(k2-1)+2) = num2str(round(med_runtimes(k1,k2)-med_runtimes(k1,length(alg_names)), n_dec), n_digits);
            end
        end
        if k2 == length(alg_names)
            res(k1,3*(k2-1)+2) = sum(runtimes(k1,:,k2)==inf);
        else
            res(k1,3*(k2-1)+3) = sum(runtimes(k1,:,k2)==inf);
        end
    end
    if k2 == length(alg_names)
        res(length(problem_names)+1,3*(k2-1)+1) = join(string(win_eq_lose), '/');
        res(length(problem_names)+1,3*(k2-1)+2) = '';
        res(length(problem_names)+2,3*(k2-1)+1) = round(mean(ranking_med_runtimes(:,k2)), n_dec);
        res(length(problem_names)+2,3*(k2-1)+2) = '';
    else
        res(length(problem_names)+1,3*(k2-1)+1) = join(string(win_eq_lose), '/');
        res(length(problem_names)+1,3*(k2-1)+2) = '';
        res(length(problem_names)+1,3*(k2-1)+3) = '';
        res(length(problem_names)+2,3*(k2-1)+1) = round(mean(ranking_med_runtimes(:,k2)), n_dec);
        res(length(problem_names)+2,3*(k2-1)+2) = '';
        res(length(problem_names)+2,3*(k2-1)+3) = '';
    end
end

%}


