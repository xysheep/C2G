function [data,label,markers] = load_mul_fcs(fdname,ctrf)
% load_mul_fcs load flow cytometry in mulfiple fcs files and return a data
% matrix as well as corresponding label vector
%   [d,l] = load_mul_fcs(fdname) load all fcs files from folder fdname into
%           d. Cell labels are assigned as 1 to number of files in the folder. 
%   [d,l] = load_mul_fcs(fdname,ctrf) load all fcs files from folder fdname
%           into d. Cell labels are assigned as 1 to number of files in the
%           folder. File name same to ctrf includes all cells whether gated
%           or ungated. Cells in ctrf that are also in any other files are 
%           removed; cells remain are ungated cell and are labeled 0.

%% Testing command: 
%% 
% Use the following scripts to test this file 
%   
%   fdname = 'test';
%   [d,l]=load_mul_fcs(fdname,'ctr.fcs');
%

%% Initiate temp variables
flist = dir(sprintf('%s/*.fcs',fdname));
flen = length(flist);
cell_label = cell(flen,1);
cell_data = cell(flen,1);
for i = 1:flen
    [td,markers] = readfcs_v2(sprintf('%s/%s',fdname,flist(i).name));
    cell_data{i} = flow_arcsinh(td,5)';
    cell_label{i} = i * ones(size(cell_data{i},1),1);
end
%% Remove duplicated cells
if exist('ctrf','var')
    flag = -1;
    for i = 1:flen
        if strcmp(ctrf,flist(i).name)
            flag = i;
            break
        end
    end
    idx = 1:flen;
    idx(idx==flag)=[];
    for i = idx
        [~,minidx] = min(pdist2(cell_data{flag},cell_data{i}));
        rmidx = minidx';
        cell_data{flag}(rmidx,:)= [];
        cell_label{flag} = zeros(size(cell_data{flag},1),1);
    end
end
%% Join the data
data = cell2mat(cell_data);
label = cell2mat(cell_label);