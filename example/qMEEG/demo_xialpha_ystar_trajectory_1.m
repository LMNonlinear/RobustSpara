%% read data
clear;clc;
input_path='Y:\qEEG\ystar\data\smooth_lg_rigel\para';
type=[];%'_MexicoChildSleep';
data_name=['smooth_lg_sigmafree_nonegalpha',type];
data_path=fullfile(input_path,data_name);

save_path = strrep(data_path, 'para', 'extract_para');
test_folder(save_path);

file = dir(fullfile(data_path,'*mat'));
% T_logs=readtable('Y:\qEEG\ystar\data\ystarlog_1563.csv');
cname={'Fp1'	'Fp2'	'F3'	'F4'	'C3'	'C4'	'P3'	'P4'	'O1'...
    'O2'	'F7'	'F8'	'T3'	'T4'	'T5'	'T6'	'Fz'	'Cz'	'Pz'};
EEG_bands={[1,4],[4,6],[6,14],[14,20]};
band_labels = {'delta','theta', 'alpha','beta'};
T_alpha = table();
T_xi= table();
T_alpha_band=table();

for isub=1:length(file)

    load(fullfile(data_path,file(isub).name));
    sub_name=xialpha.name;
    % ind = find(strcmp(T_logs.name, sub_name));
    age= xialpha(1).info.age;%unique(T_logs.age(ind));
    sex = xialpha(1).info.sex;%unique(T_logs.sex(ind));
    study = xialpha(1).info.datasetDetail;%unique(T_logs.datasetDetail(ind));
    for ichan=1:length(xialpha)
        chan_name=['ystarlog',num2str(ichan),'_',num2str(ichan)];
        r2_s=xialpha(ichan).r2_s;
        para_xi=xialpha(ichan).para_xi.kernel.para;
        Exponent=para_xi.chi;
        Offset=log10(para_xi.h);
        para_alpha=xialpha(ichan).para_alpha.kernel.para;
        CF=para_alpha.mu;
        PW=log10(para_alpha.h);
        BW=para_alpha.sigma;
        
        for iband = 1:length(EEG_bands)
            band = EEG_bands{iband};
            band_label = band_labels{iband};

            % find the index of CF in current frequency band
            idx_in_band = CF >= band(1) & CF < band(2);

            if any(idx_in_band)
                % find teh index of maxmum PW in current bands
                [max_pw, max_idx] = max(PW(idx_in_band));
                max_cf = CF(idx_in_band);
                max_bw = BW(idx_in_band);

                % save the corresponding para
                new_row_band = table({sub_name}, {study}, {chan_name}, age, {sex}, ...
                    {band_label}, max_cf(max_idx), max_pw, max_bw(max_idx), ...
                    'VariableNames', {'subName', 'study', 'chanName', 'Age', 'Sex', ...
                    'Band', 'CF', 'PW', 'BW'});
                T_alpha_band = [T_alpha_band; new_row_band];
            end
        end
        %% save result into table

        new_row_alpha = table({sub_name}, {study},{chan_name}, age, {sex}, {CF}, {PW}, {BW}, ...
            'VariableNames', {'subName','study', 'chanName', 'Age', 'Sex', 'CF', 'PW', 'BW'});
        T_alpha = [T_alpha; new_row_alpha];


        new_row_xi = table({sub_name}, {study},{chan_name}, age, {sex}, Exponent, Offset, r2_s,...
            'VariableNames', {'subName', 'study','chanName', 'Age', 'Sex', 'Exponent', 'Offset','r2'});
        T_xi = [T_xi; new_row_xi];

    end
end

table_name1=['T_alpha_hoxialpha',type];
table_name2=['T_xi_hoxialpha',type];
table_name3=['T_alpha_band_hoxialpha',type];
writetable(T_alpha, fullfile(save_path,[table_name1,'.csv']));
writetable(T_xi, fullfile(save_path,[table_name2,'.csv']));
writetable(T_alpha_band, fullfile(save_path,[table_name3,'.csv']));










% h=0.5;
% age_nature=10.^age;
% age_rep=repmat(10.^age,1,18)';
%
% figure
% plot(age_rep(:),Offset(:),'o')
%
% figure
% plot(age_rep(:),Exponent(:),'o')
%
% figure
% plot(age,Exponent(2,:)','o')
%
% xq=linspace(min(age_nature),max(age_nature),500);
% y=Offset(2,:)';
% [c,ww]= smooth(age_nature,y,'lowess');
%
% [age_sort,ind_sort]=sort(age);
% c=c(ind_sort);
% y=y(ind_sort);
% figure
% plot(age_sort,y,'o')
% hold on
% plot(age_sort,c)
%
%
%
%
%
%
%
%
%
%
%
%
%
%
% % function [yq, L, dbg] = NwSmoothInline(x, y, h, xq)
% %
% % if nargin < 4 || isempty(xq)
% %     xq = x;
% end
%
% [n, dx] = size(x);
% [nq, ~] = size(xq);
%
% if size(x, 2) > 1
%     x = reshape(x, [n, 1, dx]);
% end
%
% if size(xq, 2) > 1
%     xq = reshape(xq, [nq, 1, dx]);
% end
%
% D = x - permute(xq, [2, 1, 3]);
% Dkn = gaussian_kernel(D, h);
% yq = permute(sum(Dkn .* y, 1) ./ sum(Dkn, 1), [2, 1, 3]);
% dbg.s = sum(Dkn, 1);
%
% if nargout > 1
%     L = Dkn ./ sum(Dkn, 1).';
% end
%
% end
