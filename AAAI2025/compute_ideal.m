clear, close all


%%
pathstr = fileparts(mfilename('fullpath'));  % 本m文件所在的路径
cd(pathstr);  % 更改当前活动目录路径
addpath(genpath(fullfile(pathstr,'function')));
% filenames = getAllFiles(fullfile('Problem','Kirlik14'));
filenames = getAllFiles(fullfile('Problem','mine'));

num_workers = 5;  % 并行节点数量
if isempty(gcp('nocreate'))
    obj_pool = parpool(num_workers);
else
    obj_pool = gcp;
end

for index_file = 1 : length(filenames)
    filename = filenames{index_file};
    task_ins(index_file) = parfeval(obj_pool, @run_ins, 0, filename, 'ideal');
end


%% Monitor
flag = true;
while flag
    flag = false;
    for index_file = 1 : length(filenames)
        [~,filename_,~] = fileparts(filenames{index_file});
        switch task_ins(index_file).State
            case 'finished'
                if ~task_ins(index_file).Read
                    fetchOutputs(task_ins(index_file));
                    if isempty(task_ins(index_file).Diary)
                        disp([char(datetime) ': Task ' filename_ ' is finished (' char(task_ins(index_file).RunningDuration) ').'])
                    else
                        fprintf([char(datetime) ': Task ' filename_ ' is finished (' char(task_ins(index_file).RunningDuration) '). Note: ' task_ins(index_file).Diary])
                    end
                end
            case {'running', 'queued'}
                flag = true;
                if strcmp(task_ins(index_file).State, 'running')
                    disp([char(datetime) ': Task ' filename_ ' has been executed for ' char(task_ins(index_file).RunningDuration) '.'])
                end
            case {'failed', 'unavailable'}
                error(['Unknown error (Task ' filename_ ' is ' task_ins(index_file).State ').'])
        end
    end
    fprintf('\n')
    pause(10)
end


%% 
function run_ins(filename,algname)
    [~,filename_,~] = fileparts(filename);
    if exist(fullfile(['results_' algname],[algname '_' filename_ '.mat']), 'file')
        disp('The result already exists.');
        return
    end
    load(filename)

    options = optimoptions('intlinprog','Display','off');
    % options.ConstraintTolerance = 1e-9;
    options.ConstraintTolerance = 1e-7;
    options.IntegerTolerance = 1e-6;
    options.RelativeGapTolerance = 1e-10;
    % options.LPOptimalityTolerance = 1e-10;

    n_model_solved = zeros(p, 1);
    runtime_solver = 0;
    [y_I, n_model_solved, runtime_solver] = payoff_table(f,A,b,lb,ub,Aeq,beq,p, ...
        options, n_model_solved, runtime_solver);

    [~,~] = mkdir(['results_' algname]);
    save(fullfile(['results_' algname],[algname '_' filename_]), 'y_I', 'runtime_solver');
end


%%
function fileList = getAllFiles(directory)
% getAllFiles 递归获取指定文件夹及其所有子文件夹中的文件路径 (by GPT-4-turbo)
%
% 输入:
%   directory - 指定的起始文件夹路径
%
% 输出:
%   fileList - 所有文件的路径列表，每个文件路径作为cell数组的一个元素

    % 初始化文件列表
    fileList = {};
    
    % 获取目录中所有项的信息
    items = dir(directory);
    
    % 排除当前目录和上级目录的链接
    items = items(~ismember({items.name}, {'.', '..'}));
    
    % 遍历目录中的每一项
    for k = 1:length(items)
        % 获取当前项的完整路径
        fullPath = fullfile(items(k).folder, items(k).name);
        
        % 检查当前项是文件还是文件夹
        if items(k).isdir
            % 若是文件夹，则递归调用getAllFiles
            subFiles = getAllFiles(fullPath);
            
            % 将子文件夹中的文件路径列表合并到当前文件列表中
            fileList = [fileList; subFiles];
        else
            % 若是文件，则添加到文件列表中
            fileList = [fileList; {fullPath}];
        end
    end
end

