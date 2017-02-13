function B=modularity(A,g)
%g are index of elements whose 
k=sum(A);
m=sum(sum(A))/2;
B=A(g,g)-k(g)'*k(g)/(2*m); 