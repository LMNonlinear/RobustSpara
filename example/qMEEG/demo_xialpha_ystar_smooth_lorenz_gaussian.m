%
% path='E:\Multinational EEG Cross-Spectrum\ShareRawData';
% folders=get_folders_in_directory(path);
% files=cellfun(@(x) get_files_by_ext(fullfile(path,x),'mat'),folders,'UniformOutput',false);
% files=vertcat(files{:});
clear;clc;
close all
maxworker=20;
startmatlabpool(min(feature('numcores'),maxworker));
debug=true;
file_ystar = './data/HarMNqEEG/ystarlog_1563.csv';
tbl = readsubtable(file_ystar);
path_save='./result/qMEEG/';
task_name='smooth_lg_sigmafree_nonegalpha';
%%
subjname = unique(tbl.name);
%%
[~,idx_select]=min(unique(tbl.age,'rows','stable'));
%%
if debug
    id_subj=2;%idx_select;length(subjname);
else
    id_subj=length(subjname):-1:1;
end

%%
for isubj = id_subj%length(subjname):-1:1%idx_select% 2%%
    try
        % stbl = tbl(strcmp(tbl.name, subjname{isubj}), :);
        stbl = get_subtable(tbl, {'name', subjname{isubj}});
        stbl = sortrows(stbl, 'freq');
        header = get_spec_hearder('diag', 'ystarlog', 18);
        hos.f = get_subtable(stbl, {'name', subjname{isubj}}, 'freq').freq;
        hos.s = exp(table2array(get_subtable(stbl, {'name', subjname{isubj}}, header)));
        hos.Fs = max(hos.f) * 2;
        hos.name=subjname{isubj};
        hos.age= stbl.age(1);
        hos.sex = stbl.sex{1};
        hos.country=stbl.country{1};
        hos.device=stbl.device{1};
        hos.datasetName=stbl.datasetName{1};
        hos.datasetDetail=stbl.datasetDetail{1};
        hos.disease=stbl.disease{1};

        hos.debug=debug;
        xialpha = hoxialpha_psd(hos);
        subjid=['subj_',num2str(isubj)];
        file_save=fullfile(path_save,task_name,[subjid,'.mat']);
        test_folder(file_save);
        save(file_save,'xialpha')
        logging{isubj}='done';
        xialpha(1).show
        drawnow
        pause(0.1)
    catch
        warning(['file crash ', num2str(isubj)])
        logging{isubj}='failed';
    end
end

function xialpha = hoxialpha_psd(hos)
% delta (1–3 Hz), theta (4–7 Hz), alpha (8–12 Hz), beta (13–30 Hz), and
% gamma (30–100 Hz low and high). ->6 bands
% ks=2; kh=4;
hos = psd_smooth_nw(hos);
% idx = hos.f < 20;
% hos.f = hos.f(idx, :);
% hos.s = hos.s(idx, :);
[f, s, Fs] = deal(hos.f, hos.s, hos.Fs);
ks = 4;
kh = 8;
if hos.debug
    id_ch=1;
else
    id_ch=1:size(s,2);
end
for ich=id_ch%size(s,2):-1:1%1;1 %
    xialpha(:, ich) = BiXiAlpha(f, s(:, ich), [], ks, kh, Fs);
    xialpha(:, ich).verbose = hos.debug; %false true
    % xialpha(:, ich).svar = s(:, ich); %false true
    xialpha(:, ich).multiple_seg = false;
    xialpha(:, ich).para_xi.type = 'lorentzian';
    xialpha(:, ich).para_alpha.no_zero_peak = true;
    xialpha(:, ich).para_alpha.type = 'gaussian';
    xialpha(:, ich).regularization = 'elastic_constraint_mu_sigma_chi_xi';
    xialpha(:, ich).lambda=[1;1;1e1;1e1;1e3];%mu sigma startend noneg
    xialpha(:, ich).peak_relation='free';
    xialpha(:, ich).para_fit.StepTolerance = 1e-12;
    xialpha(:, ich).para_fit.FunctionTolerance = 1e-12;
    xialpha(:, ich).para_fit.OptimalityTolerance = 1e-12;
    xialpha(:, ich).para_fit.FiniteDifferenceStepSize = 1e-8;
    % xialpha(:, ich).para_fit.ScaleProblem = 'none';
    xialpha(:, ich).name = hos.name;
    xialpha(:, ich).info.age=hos.age;
    xialpha(:, ich).info.sex=hos.sex;
    xialpha(:, ich).info.country=hos.country;
    xialpha(:, ich).info.device=hos.device;
    xialpha(:, ich).info.datasetName=hos.datasetName;
    xialpha(:, ich).info.datasetDetail=hos.datasetDetail;
    xialpha(:, ich).info.disease=hos.disease;
end
if hos.debug
    for ich=id_ch%1:size(s,2)%1%
        tempmodel = xialpha(:, ich);
        tempmodel.fit;
        tempmodel.temp=[];
        tempmodel.model=[];
        xialpha(:, ich) = tempmodel;
    end
else
    parfor ich=id_ch%1:size(s,2)%1%
        tempmodel = xialpha(:, ich);
        tempmodel.fit;
        tempmodel.temp=[];
        tempmodel.model=[];
        xialpha(:, ich) = tempmodel;
    end
end

end

function hos = psd_smooth_nw(hos)
f_interp=linspace(min(hos.f),max(hos.f),length(hos.f)*3)';
h = 0.3; %0.4
hos.s = permute(hos.s, [1, 3, 2]);
% figure(101);clf;hold on;plot(hos.f,[hos.s(:,:,2)])
[hos.s] = NwSmoothInline(hos.f, hos.s, h,f_interp);
% [hos.s1] = LLSmooth(hos.f, hos.s, h,f_interp);
hos.f=f_interp;
% figure(101);plot(hos.f,[hos.s(:,:,2),hos.s1(:,:,2)])
hos.s = permute(hos.s, [1, 3, 2]);
end

function hos = psd_smooth_sg(hos)
% Define the polynomial order and frame size
polyOrder = 2; % Cubic polynomial
frameSize = 7; % Window size (must be odd)

% Apply the Savitzky-Golay filter
hos.s = sgolayfilt(hos.s, polyOrder, frameSize);
end

function [yq, L, dbg] = NwSmoothInline(x, y, h, xq)

if nargin < 4 || isempty(xq)
    xq = x;
end

[n, dx] = size(x);
[nq, ~] = size(xq);

if size(x, 2) > 1
    x = reshape(x, [n, 1, dx]);
end

if size(xq, 2) > 1
    xq = reshape(xq, [nq, 1, dx]);
end

D = x - permute(xq, [2, 1, 3]);
Dkn = gaussian_kernel(D, h);
yq = permute(sum(Dkn .* y, 1) ./ sum(Dkn, 1), [2, 1, 3]);
dbg.s = sum(Dkn, 1);

if nargout > 1
    L = Dkn ./ sum(Dkn, 1).';
end

end

function u = gaussian_kernel(u, b)
u = u ./ b;
u = (1 / sqrt(2 * pi)) * exp(-0.5 * sum((u .^ 2), 3));
end
