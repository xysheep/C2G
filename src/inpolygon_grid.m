function ind_in = inpolygon_grid(x,y,boundaryx,boundaryy)
if length(boundaryx) > 4
    ind_in = inpolygon(x,y,boundaryx,boundaryy);
else
    ind_in = ismember(x,boundaryx) & ismember(y,boundaryy);
end