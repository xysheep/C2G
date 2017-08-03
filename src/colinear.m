function flag = colinear(p1,p2)
if length(unique(p1)) == 1 || length(unique(p2)) == 1
    flag = 1;
elseif p1(2) ~= p1(1)
    a = (p2(2)-p2(1))./(p1(2)-p1(1));
    b = p2(1) - p1(1) *a;
    flag = sum(abs(a*p1+b - p2)) < 0.0001;
elseif p2(2) ~= p2(1)
    a = (p1(2)-p1(1))./(p2(2)-p2(1));
    b = p1(1) - p2(1) *a;
    flag = sum(abs(a*p2+b - p1)) < 0.0001;
else
    flag = 1;
end