function [STATS_STRUCT,CONFINT_STRUCT,ranef] = extract_violin_stats(rstats_table,cl_n,params)
%EXTRACT_VIOLIN_STATS Summary of this function goes here
%   Detailed explanation goes here
STR_OFFSET = [-0.1,-0.05];
STR_FONT_SIZE = 7;
CI_BAR_WIDTH = 0.15;
CI_BAR_XPOS = 0; % 0.5;
CI_BAR_LINESPECS = {'LineStyle','-','LineWidth',2,'Color','k'}; %[.9 .9 .9]
%% ===================================================================== %%
%## EXTRACT STATS INFO
if ~isfield(params,'kin_measure')
    params.kin_measure = [];
end
if isempty(params.kin_measure)
    tmp_stats = rstats_table.cluster_num==double(string(cl_n)) &...
        strcmp(rstats_table.model_char,params.model_char_int) &...
        strcmp(rstats_table.freq_band_char,params.eeg_measure) &...
        strcmp(rstats_table.group_char,params.group_char);
else
    tmp_stats = rstats_table.cluster_num==double(string(cl_n)) &...
        strcmp(rstats_table.model_char,params.model_char_int) &...
        strcmp(rstats_table.freq_band_char,params.eeg_measure) &...
        strcmp(rstats_table.kinematic_char,params.kin_measure) &...
        strcmp(rstats_table.group_char,params.group_char);
end
tmp_stats = rstats_table(tmp_stats,:);
%--
tmp_ac = strsplit(tmp_stats.anv_chars{1},',');
tmp_anv = cellfun(@(x) double(string(x)),strsplit(tmp_stats.anv_pvals{1},','));
%-- anova p-values
anvs = zeros(length(params.anv_chars_int),1);
for cc = 1:length(params.anv_chars_int)
    ind = strcmp(params.anv_chars_int{cc},tmp_ac);
    if ~isempty(ind)
        anvs(cc) = tmp_anv(ind);              
    else
        fprintf("Coefficient %s not found.\n",params.anv_chars_int{cc})
    end            
end

%## CHECK FOR SIGNIFICANT INTERACTION
if anvs(4) > 0.05
    %## USE NON-INTERACTION MODEL
    if isempty(params.kin_measure)
        tmp_stats = rstats_table.cluster_num==double(string(cl_n)) &...
            strcmp(rstats_table.model_char,params.model_char_group) &...
            strcmp(rstats_table.freq_band_char,params.eeg_measure) &...
            strcmp(rstats_table.group_char,params.group_char);
    else
        tmp_stats = rstats_table.cluster_num==double(string(cl_n)) &...
            strcmp(rstats_table.model_char,params.model_char_group) &...
            strcmp(rstats_table.freq_band_char,params.eeg_measure) &...
            strcmp(rstats_table.kinematic_char,params.kin_measure) &...
            strcmp(rstats_table.group_char,params.group_char);
    end
    tmp_stats = rstats_table(tmp_stats,:);
    %--
    tmp_ac = strsplit(tmp_stats.anv_chars{1},',');
    tmp_anv = cellfun(@(x) double(string(x)),strsplit(tmp_stats.anv_pvals{1},','));
    %-- anova p-values
    anvs = zeros(length(params.anv_chars_group),1);
    for cc = 1:length(params.anv_chars_group)
        ind = strcmp(params.anv_chars_group{cc},tmp_ac);
        if ~isempty(ind)
            anvs(cc) = tmp_anv(ind);              
        else
            fprintf("Coefficient %s not found.\n",params.anv_chars_group{cc})
        end            
    end
    tmp_cc = strsplit(tmp_stats.coeff_chars{1},',');
    tmp_fsq_chars = strsplit(tmp_stats.fsq_chars{1},',');
    tmp_ci_chars = strsplit(tmp_stats.confint_chars{1},',');
    tmp_coeffs = cellfun(@(x) double(string(x)),strsplit(tmp_stats.coeffs{1},','));
    tmp_fsq = cellfun(@(x) double(string(x)),strsplit(tmp_stats.fsq_vals{1},','));
    % tmp_em = cellfun(@(x) double(string(x)),strsplit(tmp_stats.emmeans{1},','));
    tmp_ci_lwr = cellfun(@(x) double(string(x)),strsplit(tmp_stats.confint_lwr{1},','));
    tmp_ci_upr = cellfun(@(x) double(string(x)),strsplit(tmp_stats.confint_upr{1},','));
    %-- rndm inter
    tmp_ranef_char = cellfun(@(x) string(x),strsplit(tmp_stats.ran_effs_char{1},','));
    tmp_ranef_n = cellfun(@(x) double(string(x)),strsplit(tmp_stats.ran_effs_n{1},','));
    ranef = struct('char',tmp_ranef_char, ...
        'int',tmp_ranef_n);
    %-- model coefficients       
    coeffs = zeros(length(params.coeff_chars_group),1);        
    for cc = 1:length(params.coeff_chars_group)
        ind = strcmp(params.coeff_chars_group{cc},tmp_cc);
        if ~isempty(ind)
            coeffs(cc) = tmp_coeffs(ind);              
        else
            fprintf("Coefficient %s not found.\n",params.coeff_chars_group{cc})
        end            
    end        
    %-- cohens f^2 values
    fsq_chars = strcmp(params.anv_chars_group,'(Intercept)');
    fsq_chars = params.anv_chars_group(~fsq_chars);
    fsqs = zeros(length(fsq_chars),1);
    for cc = 1:length(fsq_chars)
        ind = strcmp(fsq_chars{cc},tmp_fsq_chars);
        if ~isempty(ind)
            fsqs(cc) = tmp_fsq(ind);
        else
            fprintf("Coefficient %s not found.\n",params.anv_chars_group{cc})
        end
    end
    %-- confidence intervals
    ci_chars = strcmp(params.group_chars,'(Intercept)');
    ci_chars = params.group_chars(~ci_chars);
    cis = zeros(length(ci_chars),1,2);
    for cc = 1:length(ci_chars)
        ind = strcmp(ci_chars{cc},tmp_ci_chars);
        if ~isempty(ind)
            cis(cc,1,:) = [tmp_ci_lwr(ind),tmp_ci_upr(ind)];
        else
            fprintf("Coefficient %s not found.\n",params.group_chars{cc})
        end
    end
    %--
    % ran_effs_char = strsplit(tmp_stats.ran_effs_char{1},',');
    % ran_effs_n = cellfun(@(x) double(string(x)),strsplit(tmp_stats.ran_effs_n{1},','));
    %--
    if anvs(3) < 0.05 && anvs(3) > 0.01
        strg = '^{+}';
    elseif anvs(3) <= 0.01 && anvs(3) > 0.001
        strg = '^{++}';
    elseif anvs(3) <= 0.001
        strg = '^{+++}';
    else
        strg = '^{ns}';
    end
    %--
    if anvs(2) < 0.05 && anvs(2) > 0.01
        strs = '^{*}';
    elseif anvs(2) <= 0.01 && anvs(2) > 0.001
        strs = '^{**}';
    elseif anvs(2) <= 0.001
        strs = '^{***}';
    else
        strs = '^{ns}';
    end

    %## ASSIGN STATS
    str = {[sprintf('%sf_{s}^{2}=%1.2f    %sf_{g}^{2}=%1.2f\nR^2=%1.2f', ...
        strs,fsqs(1), ...
        strg,fsqs(2), ...
        tmp_stats.r2_c_int)],'',''};
    %--
    if anvs(3) > 0.05 && anvs(2) > 0.05
        chkd = false;
    else
        chkd = true;
    end
    CONFINT_STRUCT = struct('do_display',chkd, ...
            'y_bnds',cis, ...
            'x_vals',repmat(CI_BAR_XPOS,[3,1]), ...
            'errbar_struct',struct('line_specs',{CI_BAR_LINESPECS}, ...
                'err_bar_width',CI_BAR_WIDTH));
    % regl =[[coeffs(1),coeffs(2)]; ...
    %         [coeffs(1)+coeffs(3),coeffs(2)]; ...
    %         [coeffs(1)+coeffs(4),coeffs(2)]];
    % regl = [[coeffs(1)+coeffs(3),coeffs(2)]; ...
    %         [coeffs(1)+coeffs(4),coeffs(2)]; ...
    %         [coeffs(1),coeffs(2)]];
    % params.g_coeff_inds(:,2) = 0;
    % regl = zeros(size(params.g_coeff_inds,1),2);    
    % for i = 1:size(params.g_coeff_inds,1)
    %     %-- intercept
    %     if params.g_coeff_inds(i,1) ~= 0
    %         regl(i,1) = coeffs(1)+coeffs(params.g_coeff_inds(i,1));
    %     else
    %         regl(i,1) = coeffs(1);
    %     end
    %     %-- slope
    %     if params.g_coeff_inds(i,2) ~= 0
    %         regl(i,2) = coeffs(2)+coeffs(params.g_coeff_inds(i,2));
    %     else
    %         regl(i,2) = coeffs(2);
    %     end
    % end
    %--
    % regl = [[coeffs(1)+coeffs(3),coeffs(2)]; ... %lvl 1
    %     [coeffs(1)+coeffs(4),coeffs(2)]; ...
    %     [coeffs(1)-coeffs(3)-coeffs(4),coeffs(2)]]; %lvl 1=3
    %--
    regl = [[coeffs(1)+coeffs(3),coeffs(2)]; ... %lvl 1
        [coeffs(1)+coeffs(4),coeffs(2)]; ...
        [coeffs(1)-coeffs(3)-coeffs(4),coeffs(2)]]; %lvl 1=3
    tssc = struct('var_type','continuous', ...
        'anova',anvs(2), ...
        'multc_pvals',[],...
        'multc_pvals_pairs',[],...
        'regress_line',regl, ... 
        'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
        'regress_xvals',(0:5)*0.25,... % continuous predictors
        'order',{{}});
    tssg = struct('var_type','categorical', ...
        'anova',anvs(3), ...
        'multc_pvals',[],...
        'multc_pvals_pairs',[],...
        'regress_line',[],... 
        'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
        'regress_xvals',0,... % continuous predictors
        'order',{params.group_order});
    STATS_STRUCT = struct('sig_levels',[0.05,0.01,0.001],...
        'cond_stats',tssc, ...
        'group_stats',tssg, ...
        'stats_char',struct('str',{str}, ...
            'offsets',STR_OFFSET, ...
            'font_size',STR_FONT_SIZE, ...
            'do_display',true));
else
    %## USE INTERACTION MODEL
    tmp_cc = strsplit(tmp_stats.coeff_chars{1},',');
    tmp_fsq_chars = strsplit(tmp_stats.fsq_chars{1},',');
    tmp_ci_chars = strsplit(tmp_stats.confint_chars{1},',');
    tmp_coeffs = cellfun(@(x) double(string(x)),strsplit(tmp_stats.coeffs{1},','));
    tmp_fsq = cellfun(@(x) double(string(x)),strsplit(tmp_stats.fsq_vals{1},','));
    % tmp_em = cellfun(@(x) double(string(x)),strsplit(tmp_stats.emmeans{1},','));
    tmp_ci_lwr = cellfun(@(x) double(string(x)),strsplit(tmp_stats.confint_lwr{1},','));
    tmp_ci_upr = cellfun(@(x) double(string(x)),strsplit(tmp_stats.confint_upr{1},','));
    %-- rndm inter
    tmp_ranef_char = cellfun(@(x) string(x),strsplit(tmp_stats.ran_effs_char{1},','));
    tmp_ranef_n = cellfun(@(x) double(string(x)),strsplit(tmp_stats.ran_effs_n{1},','));
    ranef = struct('char',tmp_ranef_char, ...
        'int',tmp_ranef_n);
    %-- model coefficients
    % coeff_chars = strcmp(params.group_chars,'(Intercept)');
    % coeff_chars = params.group_chars(~coeff_chars);    
    % % params.coeff_chars_unmix = {'','group_char1','group_char2';'OHFA','OLFA','YA'};
    % params.coeff_chars_unmix = {'(Intercept)','speed_cond_num','group_char1','group_char2','speed_cond_num:group_char1','speed_cond_num:group_char2'; ...
    %     'OHFA','OHFA','OLFA','YA','OLFA','YA'};
    % % coeffs = zeros(length(params.coeff_chars_int),1);
    % coeffs = zeros(length(params.group_chars),length(params.coeff_chars_int));
    % for cc = 1:length(params.coeff_chars_int)
    %     ind = strcmp(params.coeff_chars_int{cc},tmp_cc);
    %     %-- unmix
    %     % unind = cellfun(@(x) contains(tmp_cc{ind},x) && ~isempty(x),params.coeff_chars_unmix(1,:));        
    %     % bind = cellfun(@isempty,params.coeff_chars_unmix(1,:));
    %     %--
    %     unind = cellfun(@(x) strcmp(tmp_cc{ind},x),params.coeff_chars_unmix(1,:));
    %     if any(unind) && ~isempty(ind)
    %         unindun = strcmp(params.coeff_chars_unmix{2,unind},coeff_chars);
    %         % tmp = tmp_coeffs(ind);
    %         coeffs(unindun,cc) = tmp_coeffs(ind);
    %     % elseif ~isempty(ind)
    %     %     unindun = strcmp(params.coeff_chars_unmix{2,bind},coeff_chars);
    %     %     coeffs(1,cc) = tmp_coeffs(ind); 
    %     else
    %         error("Coefficient %s not found.\n",params.coeff_chars_int{cc});            
    %     end        
    % end

    %## EXTRACT MODEL PARAMS
    %-- model coefficients
    coeffs = zeros(length(params.coeff_chars_int),1);        
    for cc = 1:length(params.coeff_chars_int)
        ind = strcmp(params.coeff_chars_int{cc},tmp_cc);
        if ~isempty(ind)
            coeffs(cc) = tmp_coeffs(ind);              
        else
            fprintf("Coefficient %s not found.\n",params.coeff_chars_group{cc})
        end            
    end  
    %-- cohens f^2 values
    fsq_chars = strcmp(params.anv_chars_int,'(Intercept)');
    fsq_chars = params.anv_chars_int(~fsq_chars);
    fsqs = zeros(length(fsq_chars),1);
    for cc = 1:length(fsq_chars)
        ind = strcmp(fsq_chars{cc},tmp_fsq_chars);
        if ~isempty(ind)
            fsqs(cc) = tmp_fsq(ind);
        else
            fprintf("Coefficient %s not found.\n",params.anv_chars_int{cc})
        end
    end
    %-- confidence intervals
    ci_chars = strcmp(params.group_chars,'(Intercept)');
    ci_chars = params.group_chars(~ci_chars);
    cis = zeros(length(ci_chars),1,2);
    for cc = 1:length(ci_chars)
        ind = strcmp(ci_chars{cc},tmp_ci_chars);
        if ~isempty(ind)
            cis(cc,1,:) = [tmp_ci_lwr(ind),tmp_ci_upr(ind)];
        else
            fprintf("Coefficient %s not found.\n",params.group_chars{cc})
        end
    end
    %--
    if anvs(4) < 0.05 && anvs(4) > 0.01
        stri = '*';
    elseif anvs(4) <= 0.01 && anvs(4) > 0.001
        stri = '**';
    elseif anvs(4) <= 0.001
        stri = '***';
    else
        stri = '^{ns}';
    end
    %## ASSIGN STATS
    str = {sprintf('%sf_{s:g}^{2}=%1.2f\nR^2=%1.2f', ...
        stri,fsqs(3),tmp_stats.r2_c_int),'',''};
    % txt_sz = 9;
    % offs = [-0.1,-0.05];
    %--
    % str = {sprintf('%sm_{ya}=%1.2f  m_{ohf}=%1.2f  m_{olf}=%1.2f\nR^2=%1.2f', ...
    %     stri,coeffs(2), ...
    %     coeffs(2)+coeffs(5), ...
    %     coeffs(2)+coeffs(6), ...
    %     tmp_stats.r2_c_int),'',''};
    % txt_sz = 9;
    % offs = [-0.11,-0.05];
    %--
    % str = {sprintf('%sy=(%1.1f)x+(%1.1f)\nR^2=%1.2f',stri,coeffs(2),coeffs(1),tmp_stats.r2_c_int), ...
    %     sprintf('y=(%1.1f)x+(%1.1f)',coeffs(2)+coeffs(5),coeffs(1)+coeffs(3)), ...
    %     sprintf('y=(%1.1f)x+(%1.1f)',coeffs(2)+coeffs(6),coeffs(1)+coeffs(4))};
    % txt_sz = 7;
    % offs = [-0.11,-0.05];
    %--
    if anvs(4) > 0.05 && anvs(3) > 0.05
        chkd = false;
    else
        chkd = true;
    end
    CONFINT_STRUCT = struct('do_display',chkd, ...
            'y_bnds',cis, ...
            'x_vals',repmat(CI_BAR_XPOS,[3,1]), ...
            'errbar_struct',struct('line_specs',{CI_BAR_LINESPECS}, ...
                'err_bar_width',CI_BAR_WIDTH));
    % regl = [[coeffs(1)+coeffs(3),coeffs(2)]; ...
    %         [coeffs(1)+coeffs(4),coeffs(2)]; ...
    %         [coeffs(1),coeffs(2)]];
    % regl = [[coeffs(1)+coeffs(3),coeffs(2)+coeffs(5)]; ...
    %         [coeffs(1)+coeffs(4),coeffs(2)+coeffs(6)]; ...
    %         [coeffs(1),coeffs(2)]];
    % regl = [[coeffs(1)+coeffs(3),coeffs(2)]; ...
    %         [coeffs(1)+coeffs(4),coeffs(2)]; ...
    %         [coeffs(1),coeffs(2)]];
    % params.g_coeff_inds(:,2) = 0;
    % regl = zeros(size(params.g_coeff_inds,1),2);    
    % for i = 1:size(params.g_coeff_inds,1)
    %     %-- intercept
    %     if params.g_coeff_inds(i,1) ~= 0
    %         regl(i,1) = coeffs(1)+coeffs(params.g_coeff_inds(i,1));
    %     else
    %         regl(i,1) = coeffs(1);
    %     end
    %     %-- slope
    %     if params.g_coeff_inds(i,2) ~= 0
    %         regl(i,2) = coeffs(2)+coeffs(params.g_coeff_inds(i,2));
    %     else
    %         regl(i,2) = coeffs(2);
    %     end
    % end
    % regl = [[coeffs(1)+coeffs(3)-coeffs(4),coeffs(2)+coeffs(5)-coeffs(6)]; ... %lvl 1
    %     [coeffs(1)+coeffs(3)+coeffs(4),coeffs(2)+coeffs(5)+coeffs(6)]; ...
    %     [coeffs(1)-coeffs(3)-coeffs(4),coeffs(2)-coeffs(5)-coeffs(6)]]; %lvl 1=3
    regl = [[coeffs(1)+coeffs(3),coeffs(2)+coeffs(5)]; ... %lvl 1
        [coeffs(1)+coeffs(4),coeffs(2)+coeffs(6)]; ...
        [coeffs(1)-coeffs(3)-coeffs(4),coeffs(2)-coeffs(5)-coeffs(6)]]; %lvl 1=3
    tssc = struct('var_type','continuous', ...
        'anova',anvs(4), ...
        'multc_pvals',[],...
        'multc_pvals_pairs',[],...
        'regress_line',regl, ...
        'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
        'regress_xvals',(0:5)*0.25,... % continuous predictors
        'order',{{}});
    tssg = struct('var_type','categorical', ...
        'anova',anvs(4), ...
        'multc_pvals',[],...
        'multc_pvals_pairs',[],...
        'regress_line',[],... 
        'line_type',{'best_fit'},... % ('best_fit' | 'means') % continuous predictors
        'regress_xvals',0,... % continuous predictors
        'order',{params.group_order});
    STATS_STRUCT = struct('sig_levels',[0.05,0.01,0.001],...
        'cond_stats',tssc, ...
        'group_stats',tssg, ...
        'stats_char',struct('str',{str}, ...
            'offsets',STR_OFFSET, ...
            'font_size',STR_FONT_SIZE, ...
            'do_display',true));
end

end

