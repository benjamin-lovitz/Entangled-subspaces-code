loadPackage "gfanInterface"

R=QQ[t_1..t_6,x_(0,0,0)..x_(2,2,1)]
c={0,1,1,1};
d={1,1,2,0};
r={1,1,1,0};
s={0,2,1,1};
S=ideal(x_(0,0,0)-t_1,x_(0,0,1)-t_2,x_(0,1,0)-t_3*c#0-t_5*d#0,x_(0,1,1),x_(0,2,0)-t_4*r#0-t_6*s#0,x_(0,2,1)-t_3*c#3-t_5*d#3,x_(1,0,0)-t_3*c#1-t_5*d#1,x_(1,0,1),x_(1,1,0)-t_1,x_(1,1,1)-t_2,x_(1,2,0),x_(1,2,1)-t_4*r#3-t_6*s#3,x_(2,0,0)-t_4*r#1-t_6*s#1,x_(2,0,1)-t_3*c#2-t_5*d#2,x_(2,1,0),x_(2,1,1)-t_4*r#2-t_6*s#2,x_(2,2,0)-t_1,x_(2,2,1)-t_2);

 -- S is a parameterization of the 6-dimensional subspace (including coefficients of each basis element, t_1..t_6)



--S'=eliminate(S,toList(t_1..t_6))
S'=gfanBuchberger(S);					   -- S' finds the reduce Grobner basis of S (which is just row reduction in this case)

s'=drop(toList S' #1,-6);				   -- s' is S' with the terms containing t_1..t_6 removed (this means our desired 6-dimensional subspace is V(s')

R=QQ[x_(0,0,0)..x_(2,2,1)];				   -- Redefine R to remove auxilliary variables t_1..t_6, and redefine S to be <s'>, an ideal of our newly defined ring R.
S=sub(ideal(s'),R);



col1=apply(flatten(table(3,2,(i,j)->x_(i,0,j)..x_(i,2,j))),toList); -- columns of 3x6 matrix formed by elements of the subspace
col2=apply(flatten(table(3,2,(i,j)->x_(0,i,j)..x_(2,i,j))),toList); -- columns of the partial transpose		

n1=transpose matrix col1;        
n2=transpose matrix col2;

D=minors(3,n1)+minors(3,n2);       -- minors of 3x6 matrix and its partial transpose
radical(S+D)			   -- radical(S+D)=(0), as expected (so S is 2-entangled)


T=QQ[t_1..t_6,x_(0,0,0,0)..x_(1,1,1,1)];
Q=ideal(x_(0,0,0,0)-t_6,x_(0,0,0,1)-t_1,x_(0,1,0,0)-t_3-t_4,x_(0,1,0,1),x_(0,0,1,0)-t_2-t_4,x_(0,0,1,1)-t_5-2*t_6,x_(0,1,1,0)-t_1,x_(0,1,1,1)-t_4,x_(1,0,0,0)-t_3,x_(1,0,0,1)-t_2,x_(1,1,0,0)-t_5-t_6,x_(1,1,0,1)-t_1-t_3,x_(1,0,1,0),x_(1,0,1,1)-t_3-t_4,x_(1,1,1,0)-t_2,x_(1,1,1,1)-t_5);
Q'=gfanBuchberger(Q);
q'=drop(toList Q' #1,-6);
T=QQ[x_(0,0,0,0)..x_(1,1,1,1)];
Q=sub(ideal(q'),T);

dol1=apply(flatten(table(2,2,(i,j)->(x_(0,i,0,j)..x_(1,i,1,j)))),toList);        -- columns of 4x4  matrix formed by elements of the subspace, and the other two flattenings
dol2=apply(flatten(table(2,2,(i,j)->(x_(0,i,j,0)..x_(1,i,j,1)))),toList);
dol3=apply(flatten(table(2,2,(i,j)->(x_(i,j,0,0)..x_(i,j,1,1)))),toList);

m1=transpose matrix dol1;
m2=transpose matrix dol2;
m3=transpose matrix dol3;

J=minors(3,m1)+minors(3,m2)+minors(3,m3);


radical(Q+J)          --radical(Q+J)=(0) as expected (so Q is 2-entangled)

