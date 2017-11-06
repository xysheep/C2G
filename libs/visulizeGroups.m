function visulizeGroups(x,y,groups,col)
%ngroup=length(unique(groups));
uni_g = unique(groups)';
for i=1:length(uni_g)
    %if i>0
    if exist('col' , 'var')
        jplot(x(groups==uni_g(i)),y(groups==uni_g(i)),'.','Color',col(i,:));
    else
        jplot(x(groups==uni_g(i)),y(groups==uni_g(i)),'.');
    end
    hold on
    
    %end
end