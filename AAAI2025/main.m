clear, close all

%% User Settings
algname = 'NPDA';
num_workers = 5;  % set properly according to your computer
path_instances = fullfile('Problem','ins_by_ins_generator');  % file path of test instances

%% Config
pathstr = fileparts(mfilename('fullpath'));
cd(pathstr);
addpath(genpath(fullfile(pathstr,'function')));
addpath(genpath(fullfile(pathstr,'MatlabSCIPInterface-master')));
filenames = getAllFiles(path_instances);

if isempty(gcp('nocreate'))
    obj_pool = parpool(num_workers);
else
    obj_pool = gcp;
end

for index_file = 1 : length(filenames)
    filename = filenames{index_file};
    task_ins(index_file) = parfeval(obj_pool, @run_ins, 0, filename, algname);
end

%% Monitor  Terminate execution using "cancel(task_ins)"
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


%% Function
function run_ins(filename,algname)
% Run an instance

    [~,filename_,~] = fileparts(filename);
    if exist(fullfile(['results_' algname],[algname '_' filename_ '.mat']), 'file')
        disp('The result already exists.');
        return
    end

    load(filename)

    tmp = strfind(filename_,'_');
    problem_name = filename_(1:(tmp(1)-1));
    switch problem_name
        case {'KP', 'AP', 'ILP', 'SPP'}
            min_finterval = 1;  % the minimum difference between different objective function values
            runtime0_k = 3600;
            options = optimoptions('intlinprog','Display','off');
            % options.ConstraintTolerance = 1e-9;
            options.ConstraintTolerance = 1e-7;
            options.IntegerTolerance = 1e-6;
            options.RelativeGapTolerance = 1e-10;
            % options.LPOptimalityTolerance = 1e-10;
            [y_N,Y_P,n_model_solved,runtime,runtime_solver,runtime_k] = feval([algname '_MOILP'], ...
                f,A,b,lb,ub,Aeq,beq,min_finterval,p, ...
                filename, ...
                runtime0_k,options);

        case 'IQP'
            min_finterval = 1;
            runtime0_k = 3600;
            options = optiset('display','off');  % 'off','iter','full'
            options.solverOpts = {'limits/time',runtime0_k; ...
                                  'numerics/feastol',1e-7; ...
                                  'numerics/epsilon',1e-9; ...
                                  'limits/gap',0; ...
                                  'limits/absgap',0};
            [y_N,Y_P,n_model_solved,runtime,runtime_solver,runtime_k] = feval([algname '_MOIQP'], ...
                H,f,A,b,lb,ub,Aeq,beq,min_finterval,p, ...
                filename, ...
                runtime0_k,options);
    end

    [~,~] = mkdir(['results_' algname]);
    save(fullfile(['results_' algname],[algname '_' filename_]), 'y_N', 'Y_P', 'n_model_solved', 'runtime', 'runtime_solver', 'runtime_k');
end


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
