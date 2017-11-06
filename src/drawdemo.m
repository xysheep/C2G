function drawdemo(i,j,col,sub_d,sub_l,over_matrix,boundary,n_markers,tmp_k,step,markernames,n_gates,entropy,flag_seperate,exclude_cells)
%isa(boundary,'cell') %&& length(boundary) > 1
    downidx = downsample(1:size(sub_d,1),1);
    plots_row = n_markers*(n_markers-1)/2;
    plots_col = 4;
    h = subplot(plots_row,plots_col,tmp_k*4+5*(step-1)-2*(step-1.5)*1);
%     if step == 2
%         p = get(h,'Position');
%         p(1) = p(1) + 0.05;
%         %p(3) = p(3) + 0.05;
%         set(h,'Position',p);
%     end

    bds = cell2mat(boundary');

    visulizeGroups(sub_d(downidx,i),sub_d(downidx,j),sub_l(downidx),col(unique(sub_l)+1,:));
    %legend(strread(num2str(unique(sub_l)'),'%s'),'Location','best');
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    xlabel(markernames(i),'FontSize', 15);
    ylabel(markernames(j),'FontSize', 15);
    axis([min([sub_d(:,i);bds(:,1)])-2 max([sub_d(:,i);bds(:,1)])+2 min([sub_d(:,j);bds(:,2)])-2 max([sub_d(:,j);bds(:,2)])+2 ])
    if tmp_k == 0
        title(sprintf('Scatter Plot of\n Cell Populations'),'FontSize', 12)
    end
    subplot(plots_row,plots_col,tmp_k*4+5*(step-1)-2*(step-1.5)*2);
    imagesc(over_matrix,[0 max(over_matrix(:))]);
    set(gca,'XTick',1:length(unique(sub_l)));
    set(gca,'XTickLabel',unique(sub_l(sub_l~=0))');
    set(gca,'YTick',1:length(unique(sub_l)));
    set(gca,'YTickLabel',unique(sub_l(sub_l~=0))');
    set(gca,'fontsize',12);%'1:length(unique(sub_l)));
    if tmp_k == 0
        title(sprintf('Compute Overlap\n between Populations'),'FontSize', 12)
    end
    subplot(plots_row,plots_col,tmp_k*4+5*(step-1)-2*(step-1.5)*3);
    m = mcl(over_matrix);
    m = round(m,3);
    clustm = zeros(size(m,1),1)';
    for mi = 1:size(m,1)
        [~,clustm(mi)] = max(m(:,mi));
    end
    displayclasify(clustm,unique(sub_l(sub_l~=0)),col);
    if tmp_k == 0
        title(sprintf('MCL Clustering\n of Populations'),'FontSize', 12)
    end
    h = subplot(plots_row,plots_col,tmp_k*4+5*(step-1)-2*(step-1.5)*4);
%     if step == 1
%         p = get(h,'Position');
%         p(1) = p(1) + 0.025;
%         %p(3) = p(3) + 0.05;
%         set(h,'Position',p);
%     end
    tmp_l = sub_l;
    visulizeGroups(sub_d(downidx,i),sub_d(downidx,j),tmp_l(downidx),col(unique(sub_l)+1,:));
    for i_bon = 1:length(boundary)
        plot(boundary{i_bon}(:,1),boundary{i_bon}(:,2),'k-','LineWidth',3);
    end
    xlabel(markernames(i),'FontSize', 15);
    ylabel(markernames(j),'FontSize', 15);
    axis([min([sub_d(:,i);bds(:,1)])-2 max([sub_d(:,i);bds(:,1)])+2 min([sub_d(:,j);bds(:,2)])-2 max([sub_d(:,j);bds(:,2)])+2 ])
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    if (n_gates>1 || (flag_seperate && exclude_cells>50) )
        str_entropy = sprintf('Entropy = %.3f',entropy);
    else
        str_entropy = 'No separation';
    end
    if tmp_k == 0
        title(sprintf('Find Best Marker Pair\n%s',str_entropy),'FontSize', 12);
    else
        title(str_entropy,'FontSize', 12);
    end
end