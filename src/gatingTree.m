classdef gatingTree < handle
%% This object is used to record gating hierarchy of flow cytometry data
    properties
        cell_idx = {[]}% index of cells in this gate
        dimpair = {[]}
        ignore_perc = {[]}
        boundary = {[]}
        main_member = {[]}% labels of cells that are true positive in this gate
        numNode = 1;
        parents    % node_id of last gate. 0 means all data
        buffersize = 10;% avoid change the buffer size in each loop
    end
%% Use the Following scripts to test the object
% 
%   
%   t = gatingTree(4);
%   t.addnode(1,1:5);
%   t.addnode(1,2:6);
%   t.addnode(3,1:5);
%   t.addnode(3,1:5);
%   t.plottree();
%
%%  Functions associated with this object
    methods
        function obj = gatingTree(num_cells,labels, buffersize)
            % Pre-assign 100 node to the tree; If Number of node exceed the
            % buffersize, double the buffersize.
            if nargin > 1 
                obj.buffersize = buffersize;
            end
            obj.cell_idx = {1:num_cells};
            obj.main_member = {labels(labels~=0)};
            obj.ignore_perc = {0};
            obj.parents = zeros(1,obj.buffersize);
        end
        function obj = addnode(obj,parent,cells,main_members,ignore_perc,boundary)
            if obj.numNode == obj.buffersize
            % If Number of nodes close to the buffersize, double the
            % buffersize
                obj.buffersize = obj.buffersize * 2;
                obj.parents = [obj.parents zeros(1,obj.numNode)];
            end
            obj.parents(obj.numNode+1) = parent;
            obj.cell_idx{obj.numNode+1} = cells;
            obj.dimpair{obj.numNode+1} = [];
            obj.main_member{obj.numNode+1} = main_members;
            obj.ignore_perc{obj.numNode+1} = ignore_perc;
            obj.boundary{obj.numNode+1} = boundary;
            obj.numNode = obj.numNode +1;
        end
        function newtree = applygates(obj,newdata)
            newtree = obj;
            newtree.cell_idx{1} = 1:size(newdata,1);
            for i = 2:obj.numNode
                parent_idx = newtree.cell_idx{obj.parents(i)};
                tmpdimpair = obj.dimpair{obj.parents(i)};
                tmpboundary = obj.boundary{i};
                newtree.cell_idx{i} = parent_idx(inpolygon(...
                    newdata(parent_idx,tmpdimpair(1)),newdata(parent_idx,tmpdimpair(2)),...
                    tmpboundary(:,1),tmpboundary(:,2)));
            end
        end
        function flag = isleaf(obj,nodeID)
            % This function check if certain node is a leaf
            flag = ~ismember(nodeID,obj.parents);
        end
        function child_IDs = children(obj,nodeID)
            % This function return node_id of all children
            child_IDs = find(obj.parents == nodeID);
        end
        function obj = setdim(obj,nodeID,pairs)
            % This function set which is the best dimension to draw gates
            obj.dimpair{nodeID} = pairs;
        end
        function plottree(obj)
            % This function plot the tree structure
            treeplot(obj.parents(1:obj.numNode));
            [x,y] = treelayout(obj.parents(1:obj.numNode));
            x=x';y=y';
            name1 = cellstr(num2str((1:obj.numNode)'));
            text(x(:,1), y(:,1), name1, ...
                'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',20)
        end
        function visulize_gating_sequence(obj,data,markers,n_lines,small_cluster,fontsize)
            if ~exist('n_lines','var')
                n_lines = 2;
            end
            if ~exist('fontsize','var')
                fontsize = 20;
            end
            if isempty(markers)
                markers =  strread(num2str(1:size(data,2)),'%s');
            end
            node_can_show = find(~cellfun(@isempty,obj.dimpair));
            % Condition some figures not shown:
            % 1. Only one gate drawn on this step
            % 2. The only drawn gate only exclude small clusters (defined
            % by variable "small_clusters"
            node_to_show = node_can_show;
            if exist('small_cluster','var')
                for n_i = length(node_to_show):-1:1
                    if length(find(obj.parents==node_to_show(n_i)))==1
                        self_i = node_to_show(n_i);
                        child_i = find(obj.parents==node_to_show(n_i),1);
                        n_exclude = length(obj.cell_idx{self_i})-length(obj.cell_idx{child_i});
                        if n_exclude < small_cluster
                            node_to_show(n_i) = [];
                        end
                    end
                end
            end
            figure;
            set(gcf,'Position',[100,100,320*ceil(length(node_to_show)/n_lines),320*n_lines])
            flag = 0;
            for n_i = 1:length(node_to_show)
                n_id = node_to_show(n_i);
                if n_lines> 7 && n_i > length(node_to_show)/2
                    if flag == 0
                        figure;
                        set(gcf,'Position',[100,100,1200,600])
                        flag = 1;
                    end
                    subplot(ceil(n_lines/2),ceil(length(node_to_show)/n_lines),n_i - floor(length(node_to_show)/2));
                elseif n_lines> 7
                    subplot(ceil(n_lines/2),ceil(length(node_to_show)/n_lines),n_i);
                else
                    subplot(n_lines,ceil(length(node_to_show)/n_lines),n_i);
                end
                title(sprintf('Node %d',n_id),'FontSize',fontsize);
                sub_d = data(obj.cell_idx{n_id},:);
                %disp(obj.dimpair{n_id})
                i = obj.dimpair{n_id}(1);j = obj.dimpair{n_id}(2);
                scatplot(sub_d(:,i),sub_d(:,j),'voronoi',[],100,5,1,4);
                xlabel(markers{i},'FontSize',fontsize);
                ylabel(markers{j},'FontSize',fontsize);
                hold on
                gate_mem = find(obj.parents == n_id);
                for g_id = gate_mem
                    %disp(g_id)
                    if isempty(obj.boundary{g_id})% || isequal(obj.boundary{g_id},0)
                        xg = data(obj.cell_idx{g_id},i);
                        yg = data(obj.cell_idx{g_id},j);
                        [boundaryx,boundaryy] = findboundary(xg,yg); 
                    else
                        boundaryx = obj.boundary{g_id}(:,1);
                        boundaryy = obj.boundary{g_id}(:,2);
                    end
                    plot(boundaryx,boundaryy,'r','LineWidth',2);
                    t=text(mean(boundaryx),mean(boundaryy),sprintf('Node %d',g_id),'HorizontalAlignment','center');
                    t.FontSize = fontsize;
                    t.FontWeight = 'bold';
                end
            end
            %subplot(n_lines,ceil((1+length(node_to_show))/n_lines),length(node_to_show)+1);
            figure('Position',[100,100,320,320])
            obj.plottree()
        end
        function [outtable,nmi] = show_f_score(obj,label)
            idx = 1:obj.numNode;
            leaf_idx = idx(~ismember(idx,unique(obj.parents)) & cellfun(@length,obj.cell_idx)>0);
            outtable = zeros(length(unique(label))-1,6);
            k = 0;
            for i = 1:length(leaf_idx)
                tar_pop = obj.main_member{leaf_idx(i)};
                tar_pop = tar_pop(ismember(tar_pop,unique(label)));
                for j = 1:length(tar_pop)
                    k=k+1;
                    outtable(k,1) = tar_pop(j);
                    outtable(k,2) = leaf_idx(i);
                    outtable(k,3) = sum(label(obj.cell_idx{leaf_idx(i)}) == tar_pop(j));%True positive
                    outtable(k,4) = sum(label(obj.cell_idx{leaf_idx(i)}) ~= tar_pop(j));%False positive
                    outtable(k,5) = sum(~ismember(find(label==tar_pop(j)),obj.cell_idx{leaf_idx(i)}));%False negative
                end
            end
            outtable(:,6) = 2*outtable(:,3)./(2*outtable(:,3)+outtable(:,4)+outtable(:,5));
            [~,pop_order] = sort(outtable(:,1));
            outtable = outtable(pop_order,:);
            fprintf('Population\tGate\t    TP\t    FP\t    FN\tF-score\n');
            fprintf('%10d\t%4d\t%6d\t%6d\t%6d\t  %.3f\n',outtable');
            nmi = nmi_gate(label,obj.cell_idx(leaf_idx));
            fprintf('Normalized Mutual Information = %.3f\n',nmi);
        end
        function can_merge(obj)
            idx = 1:obj.numNode;
            leaf_idx = idx(~ismember(idx,unique(obj.parents)));
            adj = zeros(length(leaf_idx));
            adj(1:(length(leaf_idx)+1):end) = cellfun(@length,obj.cell_idx(leaf_idx));
            for i = 1 : length(leaf_idx) - 1
                for j = i+1 : length(leaf_idx)
                    adj(i,j) = length(intersect(obj.cell_idx{leaf_idx(i)},obj.cell_idx{leaf_idx(j)}));
                    adj(i,j) = adj(i,j) .^2 ./ (adj(i,i)*adj(j,j));
                    adj(j,i) = adj(i,j);
                end
            end
            adj(1:(length(leaf_idx)+1):end) = 1;
            m_c = mcl(adj);
            [~,mx] = max(m_c);
            for i = unique(mx)
                merge_pop = cell2mat(obj.main_member(leaf_idx(mx==i)));
                if length(merge_pop)> 1 && obj.check_parents(leaf_idx(mx==i))
                    fprintf('Population %s can be merged, ',sprintf('%d ' ,merge_pop(:)));
                    fprintf('they are in gate %s.\n',sprintf('%d ' ,leaf_idx(mx==i)));
                end
            end
            %figure;imagesc(adj)
        end
        function can_merge = check_parents(obj,node_idx)
            parent_idx = unique(node_idx);
            can_merge = 0;
            while length(unique(parent_idx)) > 1%sum(obj.parents == parent_idx(i)) < 2 
                %lower_idx = (parent_idx~=min(parent_idx));
                %parent_idx(lower_idx) = obj.parents(parent_idx(lower_idx))
                tmp_parent_idx = parent_idx;
                for p = unique(obj.parents(parent_idx))
                    if sum(obj.parents(parent_idx) == p) > 1 % multiple node have save parents, check if they are only children
                        if sum(obj.parents(parent_idx) == p) ~= sum(obj.parents == p) 
                            return
                        else
                            parent_idx = [p parent_idx(obj.parents(parent_idx) ~= p)];
                        end
                    elseif sum(obj.parents == p) == 1 % only child, just climb up
                        parent_idx(obj.parents(parent_idx) == p) = obj.parents(parent_idx(obj.parents(parent_idx) == p));
                    end
                end
                if isequal(tmp_parent_idx,parent_idx)
                    return
                end
            end
            if parent_idx(1) ~= 1
                can_merge = 1;
            end
            %can_merge = (length(unique(parent_idx)) == 1 && sum(obj.parents == parent_idx(1)));
        end
    end
    
end
            