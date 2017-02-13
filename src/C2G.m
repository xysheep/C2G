function T = C2G(d,l,ori_l,density,varargin)%ig_ratio,markernames,col)
% C2G perform the analysis return a gatingTree object that store the
% obtained gating hierarchy.
%       T = C2G(d,l,ori_l,...) "d" is the M-by-N data matrix where M is the
%       number of markers. "l" and "ori_l" are M-by-1 matrix represent cell
%       labels after and before pre-cluster. 
%       T = C2G(d,l,ori_l,density,...) "density" is the precomputed local
%       density. It can be computed outside C2G use function
%       "compute_density". If this parameter is not provided, C2G will
%       compute local density later.
% Optional parameter:
%       ignore_ratio: percentage of low density cells ignored when compute
%       overlap between different populations. If the data is expected to
%       be highly noisy, increase this value for better performance.
%       Default is 0.05. 

n_markers = size(d,2);
% Compute local density
if ~exist('density','var')
    fprintf('Local density not provided, computing...\n');tic;
    density = compute_density(d,l);toc;
end

% Initiate other parameters
pnames = { 'ignore_ratio','markernames','color'};
dflts  = { 0.05          ,[]           ,[]};
[ig_ratio,markernames,col] = internal.stats.parseArgs(pnames,dflts,varargin{:});


queue = CQueue();
queue.push(1);
T = gatingTree(length(l),unique(l),4);
step = 2;
% T is the object representing the tree of gating hierachy
while ~queue.isempty()
    node_id = queue.pop();
    cells_idx = T.cell_idx{node_id};
    %if length(unique(ori_l(ismember(1:length(l),cells_idx)' & l~=0)))>1% || gate_fdr(cells_idx,l)
    %if length(T.main_member{node_id})>1
        %tabulate(categorical(l(cells_idx)));
        % Initiate temp variables
        best_n_gates = 1;
        best_entropy = inf;
        best_pair = 1;
        best_gatelabels = cell(1);
        best_mainmem = cell(1);
        best_igperc = cell(1);
        best_boundary = cell(1);
        %best_exclude_gates =0;
        sub_d = d(cells_idx,:);
        sub_l = l(cells_idx);
        %sub_l(~ismember(sub_l,T.main_member{node_id})) = 0;
        sub_ori_l = ori_l(cells_idx);
        sub_ori_l(~ismember(sub_ori_l,T.main_member{node_id})) = 0;
        %fprintf('Perfect Entropy is %.2f\n',perfect_entropy);
        % Reason for above command: If one gate is drawn to include a
        % certain target population,this gate can still include some cells
        % from other population. To draw next gate based this, no need to
        % draw a separate gate for other populations. 
        if step>=1 && ~isempty(markernames) && ~isempty(col)
            f = figure;
            set(f,'Position',[100 100 1400 900])
        end
        
        tmp_k = 0;
        for i = 1:n_markers - 1
            for j = i+1:n_markers
                %tic
                %fprintf('i=%2d\tj=%2d\t',i,j);
                
                
                [gatelabels,main_members,ignore_perc,boundary,flag_seperate,over_matrix] = ...
                    new_bestgate(sub_d(:,i),sub_d(:,j),...
                    sub_l,unique(sub_ori_l),T.main_member{node_id},...
                    density(cells_idx,i,j),ig_ratio);
                
                entropy = new_entropy_gate(sub_ori_l,gatelabels);
                %fprintf('used_clucster = %d\tentropy=%.2f\n',length(gatelabels),entropy);
                %toc
                n_gates = length(gatelabels);exclude_cells=0;
                if n_gates ==1
                    %current_pop = length(T.main_member{node_id});
                    %gated_pop = length(main_members{1}');
                    % Termination condition is too strong here!
                    current_pop = length(sub_l);
                    gated_pop = length(sub_l(gatelabels{1}));
                    exclude_cells = current_pop - gated_pop;
                end
                if (n_gates>1 || (flag_seperate && exclude_cells>50) ) && (entropy < best_entropy ||...
                        (entropy==best_entropy && n_gates > best_n_gates))
                    best_entropy = entropy;
                    best_n_gates = n_gates;
                    %best_exclude_gates = exclude_gates;
                    best_pair = [i j];
                    best_gatelabels = gatelabels;
                    best_mainmem = main_members;
                    best_igperc = ignore_perc;
                    best_boundary = boundary;

                end
                
                if step >= 1 && ~isempty(markernames)&& ~isempty(col)%isa(boundary,'cell') %&& length(boundary) > 1
                    h = subplot(n_markers*(n_markers-1)/2,4,tmp_k*4+5*(step-1)-2*(step-1.5)*1);
                    if step == 2
                        p = get(h,'Position');
                        p(1) = p(1) + 0.05;
                        %p(3) = p(3) + 0.05;
                        set(h,'Position',p);
                    end
                    
                    
                    
                    visulizeGroups(sub_d(:,i),sub_d(:,j),sub_l,col(unique(sub_l)+1,:));
                    legend(strread(num2str(unique(sub_l)'),'%s'),'Location','best');
                    xlabel(markernames(i),'FontSize', 15);
                    ylabel(markernames(j),'FontSize', 15);
                    if tmp_k == 0
                        title(sprintf('Scatter Plot of\n Cell Populations'),'FontSize', 16)
                    end
                    subplot(n_markers*(n_markers-1)/2,4,tmp_k*4+5*(step-1)-2*(step-1.5)*2);
                    imagesc(over_matrix,[0 max(over_matrix(:))]);
                    set(gca,'XTick',1:length(unique(sub_l)));
                    set(gca,'XTickLabel',unique(sub_l(sub_l~=0))');
                    set(gca,'YTick',1:length(unique(sub_l)));
                    set(gca,'YTickLabel',unique(sub_l(sub_l~=0))');
                    if tmp_k == 0
                        title(sprintf('Compute Overlap\n between Populations'),'FontSize', 16)
                    end
                    subplot(n_markers*(n_markers-1)/2,4,tmp_k*4+5*(step-1)-2*(step-1.5)*3);
                    m = mcl(over_matrix);
                    m = round(m,3);
                    clustm = zeros(size(m,1),1)';
                    for mi = 1:size(m,1)
                        [~,clustm(mi)] = max(m(:,mi));
                    end
                    displayclasify(clustm,unique(sub_l(sub_l~=0)),col);
                    if tmp_k == 0
                        title(sprintf('MCL Clustering\n of Populations'),'FontSize', 16)
                    end
                    h = subplot(n_markers*(n_markers-1)/2,4,tmp_k*4+5*(step-1)-2*(step-1.5)*4);
                    if step == 1
                        p = get(h,'Position');
                        p(1) = p(1) + 0.025;
                        %p(3) = p(3) + 0.05;
                        set(h,'Position',p);
                    end
                    tmp_l = sub_l;
                    visulizeGroups(sub_d(:,i),sub_d(:,j),tmp_l,col(unique(sub_l)+1,:));
                    for i_bon = 1:length(boundary)
                        plot(boundary{i_bon}(:,1),boundary{i_bon}(:,2),'b-o');
                    end
                    xlabel(markernames(i),'FontSize', 15);
                    ylabel(markernames(j),'FontSize', 15);
                    if (n_gates>1 || (flag_seperate && exclude_cells>50) )
                        str_entropy = sprintf('Entropy = %.3f',entropy);
                    else
                        str_entropy = 'No separation';
                    end
                    if tmp_k == 0
                        title(sprintf('Find Best Marker Pair\n%s',str_entropy),'FontSize', 16);
                    else
                        title(str_entropy,'FontSize', 16);
                    end
                end
                tmp_k = tmp_k + 1;
            end
        end
        step = step - 1;
%         best_gatelabels
% %         fprintf('i=%2d\tj=%2d\tn=%2d\tn=%2d\t%.3f\n',best_pair,best_n_gates,best_exclude_gates,best_entropy)
%         if length(best_pair)>1
%             figure
%             disp(best_pair)
%             i = best_pair(1);j=best_pair(2);
%             scatplot(sub_d(:,i),sub_d(:,j),'voronoi',[],100,5,2,4);
%             hold on;
%             xlabel(i);
%             ylabel(j);
%             for t = 1:length(best_gatelabels)
%                 xg = sub_d(best_gatelabels{t},i);
%                 yg = sub_d(best_gatelabels{t},j);
%                 [boundaryx,boundaryy] = findboundary(xg,yg);  
%                 plot(boundaryx,boundaryy,'LineWidth',2);
%             end
%             drawnow;
%         end
        
        if length(best_pair)>1
            %fprintf('[Selected] i = %d\tj=%d\n',best_pair(1),best_pair(2));
            T.setdim(node_id,best_pair);
            for i_gate=1:length(best_gatelabels)
                T.addnode(node_id,cells_idx(best_gatelabels{i_gate}),...
                    best_mainmem{i_gate},best_igperc{i_gate},best_boundary{i_gate});
                queue.push(T.numNode);
            end
        end
   % end
end


    