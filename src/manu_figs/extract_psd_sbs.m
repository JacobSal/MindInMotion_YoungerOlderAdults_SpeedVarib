function [psd_out] = extract_psd_sbs(dat_out_struct,dat_calcs,cl_i,conds_out,groups_out)
%LOCAL_PSD_PLOT Summary of this function goes here
%   Detailed explanation goes here
%## EXTRACT PSD DATA =================================================== %%
%-- temp params
psd_proc_char = [];
%## EXTRACT DATA
switch dat_calcs{1}
    case 'cov'
        psd_dat_out1 = dat_out_struct.psd_dat;      
        psd_dat_out2 = dat_out_struct.psd_std_dat;
        psd_dat_out = psd_dat_out2./psd_dat_out1;
        psd_proc_char = [psd_proc_char,'cov'];
    case 'mean'
        psd_dat_out = dat_out_struct.psd_dat;
        psd_proc_char = [psd_proc_char,'mean'];
    case 'std'
        psd_dat_out = dat_out_struct.psd_std_dat;
        psd_proc_char = [psd_proc_char,'std'];
end

%## EXTRACT SPECIFICS OF DATA
cond_dat_out = dat_out_struct.cond_dat;
subj_dat_out = dat_out_struct.subj_dat;
group_dat_out = dat_out_struct.group_dat;

%## PSD DATA
tmp_dat = squeeze(psd_dat_out(:,:,:,cl_i)); %[subject, frequency, epoch/splice, channel/component]
tmp_dat = reshape(permute(tmp_dat,[3,1,2]),size(tmp_dat,1)*size(tmp_dat,3),size(tmp_dat,2)); %[subject x epoch/splice, frequency];
chk = all(~isnan(tmp_dat),2);
tmp_dat = tmp_dat(chk,:);
chk = ~all(tmp_dat==0,2);
tmp_dat = tmp_dat(chk,:);
% if all(chk)
%     chk = ~all(tmp_dat==0,2);
% end
% tmp_dat = tmp_dat(chk,:);
% sum(chk)

%## CONDITION DATA
tmp_cond = squeeze(cond_dat_out(:,:,cl_i)); %[subject, frequency, epoch/splice, channel/component]    
tmp_cond = reshape(permute(tmp_cond,[2,1]),[size(cond_dat_out,1)*size(cond_dat_out,2),1]);
chk = cellfun(@isempty,tmp_cond);
tmp_cond = tmp_cond(~chk);
%--
tmp_subj = squeeze(subj_dat_out(:,:,cl_i)); %[subject, frequency, epoch/splice, channel/component]    
tmp_subj = reshape(permute(tmp_subj,[2,1]),[size(subj_dat_out,1)*size(subj_dat_out,2),1]);
chk = all(~isnan(tmp_subj),2);
tmp_subj = tmp_subj(chk,:);
chk = ~all(tmp_subj==0,2);
tmp_subj = tmp_subj(chk,:);
% if all(chk)
%     chk = ~all(tmp_subj==0,2);
% end
% tmp_subj = tmp_subj(chk,:);
subjs = unique(tmp_subj);

%## GROUPING DATA
tmp_group = squeeze(group_dat_out(:,:,cl_i)); %[subject, frequency, epoch/splice, channel/component]    
tmp_group = reshape(permute(tmp_group,[2,1]),[size(group_dat_out,1)*size(group_dat_out,2),1]);
chk = cellfun(@isempty,tmp_group);
tmp_group = tmp_group(~chk);

%## GET GOOD INDICES
% conds = unique(tmp_cond);        
% group_chars = unique(tmp_group);
% indsc = cellfun(@(x) any(strcmp(x,conds_out)),conds);
% conds = conds(indsc);
% indsc = cellfun(@(x) any(strcmp(x,conds)),tmp_cond);
%-- get condition groups requested
indsc = cellfun(@(x) any(strcmp(x,conds_out)),tmp_cond);
indsg = cellfun(@(x) any(strcmp(x,groups_out)),tmp_group);
indscg = indsg & indsc;
%-- index
tmp_cond = tmp_cond(indscg,:);
tmp_subj = tmp_subj(indscg,:);
tmp_dat = tmp_dat(indscg,:);
tmp_group = tmp_group(indscg,:);

%## FORMAT DATA
psd_out = cell(length(conds_out),1);
for g_ii = 1:length(groups_out)
    indsg = strcmp(tmp_group,groups_out{g_ii});
    for c_ii = 1:length(conds_out)
        indc = strcmp(tmp_cond,conds_out{c_ii});
        tmp = nan(size(tmp_dat,2),length(subjs));
        for s_i = 1:length(subjs)
            inds = tmp_subj == subjs(s_i);
            chk = indc & inds & indsg;
            switch dat_calcs{2}
                case 'mean'
                    %-- mean
                    tmp(:,s_i) = mean(tmp_dat(chk,:),1);
                    if g_ii == 1 && c_ii == 1 && s_i == 1
                        psd_proc_char = [psd_proc_char,'mean'];
                    end
                case 'std'
                    %-- std
                    tmp(:,s_i) = std(tmp_dat(chk,:),[],1);
                    if g_ii == 1 && c_ii == 1 && s_i == 1
                        psd_proc_char = [psd_proc_char,'std'];
                    end
                case 'prctile'
                    %-- prctile
                    tmp(:,s_i) = prctile(tmp_dat(chk,:),75,1) - prctile(tmp_dat(chk,:),25,1);
                    if g_ii == 1 && c_ii == 1 && s_i == 1
                        psd_proc_char = [psd_proc_char,'prct'];
                    end
            end
        end
        % tmp = tmp(:,all(tmp ~= 0,1));
        tmp = tmp(:,all(~isnan(tmp),1));
        psd_out{c_ii,g_ii} = tmp; %tmp_dat(inds,:);
    end
end
fprintf('done. extract_psd_sbs.m\n')
end

