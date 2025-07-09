function [ mask, thres, inds, fdr_val ] = local_fdr( p_list, alpha, corrected )
% Computes the False Discovery Rate according to Benjamini and Hochberg (1995). 
% 
% Inputs: 
% p_list - list of p values
% alpha - the desired alpha threshold. Default: 0.05
% corrected - set to true if correction for dependencies is to be applied, according to Benjamini
% and Yekutieli (2001) (this is probably not the common case). 
%
% outputs:
% ind - the indexes of significant p-values within p_list
% thres - the p-value which served as the actual threshold in this test. 
% 
% Written by Edden Gerber, lab of Leon Y. Deouell, 2012
% Please send bug reports and requsts to edden.gerber@gmail.com
%

%## =================================================================== %%#
sz = size(p_list);
p_list = reshape(p_list,numel(p_list),1);
n_vals = length(p_list);
num_tests = n_vals; 
%(06/05/2025) Anonymos (JS), there was some reason that in some cases you may want to set this to
% a lower value, but I don't remember what it was. 

if nargin < 2
    alpha = 0.05;
end

if nargin < 3
    corrected = false;
end

%-- sort in descending order
p_sorted = sort(p_list,'descend');
% osig_mask = p_sorted <= alpha;

%-- generate vector of Benjamini-Hochberg critical values: BH_crit = (i/m)*q
if corrected
    comp = (num_tests:-1:1)/num_tests * alpha / sum((1:num_tests)/num_tests);
else
    comp = (num_tests:-1:1)/num_tests * alpha;
end

%-- no idea what this does
lnc = length(comp);
comp = comp((lnc-n_vals+1):lnc);

%-- find the largest ith value of BH_crits
% tmp = p_sorted' <= comp';
% ind = zeros(length(p_sorted),1);
% for i = 1:length(p_sorted)
%     ttmp = find(tmp(i,:),1,'first');
%     if ~isempty(ttmp)
%         ind(i) = ttmp;
%     else
%         ind(i) = 0;
%     end
% end
% i = max(ind);

tmp = p_sorted <= comp;
i = find(tmp,1,'first');

if isempty(i) || i == 0
    thres = 0;
else
    thres = p_sorted(i);
end

%-- reject values less than the new critical p-value
mask = (p_list<=thres);
inds = find(mask);

%-- calculate the FDR value
if sum(mask) > 0
    % fdr_val = sum((~mask).*osig_mask)/sum(mask);
    fdr_val = (thres*sum(mask))/sum(mask);
else
    fdr_val = 1;
end

%-- reshape
mask = reshape(mask,sz);

end

