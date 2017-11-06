classdef gatingTree < handle
%% This object is used to record gating hierarchy of flow cytometry data
    properties
        cell_label = {[]}% index of cells in this gate
        cell_idx = {[]}
        boundary = {[]}
        dimpair = {[]}
        numNode = 1;
        parents    % node_id of last gate. 0 means all data
        buffersize = 10;% avoid change the buffer size in each loop
        depth = {[]}
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
        function obj = gatingTree(all_label, numcell, buffersize)
            % Pre-assign 100 node to the tree; If Number of node exceed the
            % buffersize, double the buffersize.
            if nargin > 2 
                obj.buffersize = buffersize;
            end
            if nargin > 1
                obj.cell_idx = {1:numcell};
            end
            obj.cell_label{1} = all_label(all_label~=0);
            obj.depth{1} = 0;
            obj.parents = zeros(1,obj.buffersize);
        end
        function obj = addnode(obj, parent, cell_labels)
            if obj.numNode == obj.buffersize
            % If Number of nodes close to the buffersize, double the
            % buffersize
                obj.buffersize = obj.buffersize * 2;
                obj.parents = [obj.parents zeros(1,obj.numNode)];
            end
            obj.parents(obj.numNode+1) = parent;
            obj.depth{obj.numNode+1} = obj.depth{parent} + 1;
            obj.cell_label{obj.numNode+1} = cell_labels;
            obj.numNode = obj.numNode +1;
        end
        function obj = addnode1(obj,parent, cell_idx, cell_labels, boundary)
            if obj.numNode == obj.buffersize
            % If Number of nodes close to the buffersize, double the
            % buffersize
                obj.buffersize = obj.buffersize * 2;
                obj.parents = [obj.parents zeros(1,obj.numNode)];
            end
            obj.parents(obj.numNode+1) = parent;
            obj.depth{obj.numNode+1} = obj.depth{parent} + 1;
            obj.cell_idx{obj.numNode+1} = cell_idx;
            obj.cell_label{obj.numNode+1} = cell_labels;
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
        function draw_gates(obj,data,label,density)
            for i_node = 2:obj.numNode
                pid = obj.parents(i_node);
                idx = obj.cell_idx{pid};
                x = data(idx,obj.dimpair{pid}(1));
                y = data(idx,obj.dimpair{pid}(2));
                dens2d = density(:,obj.dimpair{pid}(1),obj.dimpair{pid}(2));
                subl = label(idx,:);
                unil = unique(obj.cell_label{i_node});unil(unil==0) = [];
                boundariesx = [];
                boundariesy = [];
                idx_cell = cell(length(unil),1);
                for i = 1:length(unil)
                    xg = x(subl == unil(i));
                    yg = y(subl == unil(i));
                    %fprintf('xg: %d\tyg %d\n',length(xg),length(yg));
                    %[boundaryx,boundaryy] = findboundary(xg,yg);
%                     figure
%                     hold on
%                     plot(boundaryx,boundaryy)
                    [boundaryx,boundaryy] = findboundary_outliers(xg,yg,dens2d(subl == unil(i)),0.2);
%                     plot(boundaryx,boundaryy)
%                     plot(mean(xg),mean(yg),'o','markersize',20);
%                     plotdensity(xg,yg,dens2d(subl == unil(i)))
%                     title(sprintf('selected %d\t origin %d\n',length(select_idx), length(xg)));
                    if length(boundaryx) > 2 % At least 3 points to be a polygon
                        in_idx = inpolygon(x,y,boundaryx,boundaryy);
                        fprintf('Node: %d\t Cluster: %d\t Cells: %d\n',i_node,unil(i),sum(in_idx));
                        idx_cell{i} = reshape(idx(in_idx), numel(idx(in_idx)),1);
                        if isempty(boundariesy)
                            boundariesx = boundaryx;
                            boundariesy = boundaryy;
                        else
                            [boundariesx,boundariesy] = polybool('union',boundariesx,boundariesy,boundaryx,boundaryy);
                        end
                    end                   
                end
                obj.cell_idx{i_node} = unique(cell2mat(idx_cell));
                obj.boundary{i_node} = [boundariesx,boundariesy]; 
            end
        end
        function view_gates(obj,data,markers,varargin)
            pnames = {'fontsize','n_lines','ignore_small','onepanel'};
            dflts = {20, 3, 100, false};
            [fontsize, n_lines, ignore_small, onepanel] = internal.stats.parseArgs(pnames,dflts,varargin{:});
            if isempty(markers)
                markers =  strread(num2str(1:size(data,2)),'%s');
            end
            node_to_show = find(~cellfun(@isempty,obj.dimpair));
            num_gates = histc(obj.parents,node_to_show);
            idx = false(size(node_to_show));
            idx( num_gates > 1) = true;
            parentNum = cellfun(@length,obj.cell_idx(node_to_show));
            childNum = cellfun(@length,...
                obj.cell_idx(arrayfun(@(x)find(obj.parents == x,1),node_to_show)));
            idx( num_gates == 1 & (parentNum > ignore_small + childNum)) = true;
            node_to_show = node_to_show(idx);
            figure;
            set(gcf,'Position',[50,50,320*ceil((length(node_to_show)+onepanel)/n_lines),320*n_lines])
            for n_i = 1:length(node_to_show)
                subplot(n_lines,ceil((length(node_to_show)+onepanel)/n_lines),n_i);
                n_id = node_to_show(n_i);
                title(sprintf('Node %d',n_id),'FontSize',fontsize);
                sub_d = data(obj.cell_idx{n_id},:);
                %disp(obj.dimpair{n_id})
                i = obj.dimpair{n_id}(1);j = obj.dimpair{n_id}(2);
                scatplot(sub_d(:,i),sub_d(:,j),'gd',[],100,5,1,4);
                xlabel(markers{i},'FontSize',fontsize);
                ylabel(markers{j},'FontSize',fontsize);
                hold on
                gate_mem = find(obj.parents == n_id);
                for g_id = gate_mem
                    %disp(g_id)
                    boundaryx = obj.boundary{g_id}(:,1);
                    boundaryy = obj.boundary{g_id}(:,2);
                    plot(boundaryx,boundaryy,'r','LineWidth',2);
                    t=text(mean(boundaryx),mean(boundaryy),sprintf('Node %d',g_id),'HorizontalAlignment','center');
                    t.FontSize = fontsize;
                    t.FontWeight = 'bold';
                end
            end
            %subplot(n_lines,ceil((1+length(node_to_show))/n_lines),length(node_to_show)+1);
            if onepanel
                subplot(n_lines,ceil((length(node_to_show)+onepanel)/n_lines),n_i+1);
            else
                figure('Position',[100,100,640,640])
            end
            obj.plottree()
        end
        function [outtable,nmi] = show_f_score(obj,label)
            idx = 1:obj.numNode;
            leaf_idx = idx(~ismember(idx,unique(obj.parents)) & cellfun(@length,obj.cell_idx)>0);
            outtable = zeros(length(unique(label))-1,7);
            k = 0;
            for i = 1:length(leaf_idx)
                tar_pop = obj.cell_label{leaf_idx(i)};
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
            outtable(:,7) = 100*(outtable(:,3) + outtable(:,5))./(length(label));
            [~,pop_order] = sort(outtable(:,1));
            outtable = outtable(pop_order,:);
            fprintf('Population\tGate\t    TP\t    FP\t    FN\tF-score\tPercentage\n');
            fprintf('%10d\t%4d\t%6d\t%6d\t%6d\t  %.3f\t%2.1f\n',outtable');
            nmi = nmi_gate(label,obj.cell_idx(leaf_idx));
            fprintf('Normalized Mutual Information = %.3f\n',nmi);
            fprintf('Average F-score = %.3f\n',mean(outtable(:,6)));
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
                merge_pop = cell2mat(obj.cell_label(leaf_idx(mx==i)));
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
            