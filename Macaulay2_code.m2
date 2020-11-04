-- This Macaulay2 code is supplementary material for the paper "Entangled subspaces and generic local state discrimination with
-- pre-shared entanglement" by Benjamin Lovitz and Nathaniel Johnston, available as an arXiv preprint here:  https://arxiv.org/abs/2010.02876.

-- We first verify that the 6-dimensional subspace of C^3 otimes C^3 otimes C^2 presented in Equation (12) is 2-entangled, and then do the same for the
-- 6-dimensional subspace of C^2 otimes C^2 otimes C^2 otimes C^2 presented in Equation (13).


loadPackage "gfanInterface"


-- We first verify that the 6-dimensional subspace of C^3 otimes C^3 otimes C^2 presented in Equation (12) is 2-entangled. As mentioned in the text,
-- it suffices to check that the flattening rank of every element of the subspace is at least three. Since the coordinates of the basis elements of the
-- subspace are integers, and the determinantal equations are also integers, it suffices to work over the field of rational numbers (we expand on this
-- point inside the code).

R=QQ[t_1..t_6,x_(0,0,0)..x_(2,2,1)];
c={0,1,1,1}; -- delta in the text
d={1,1,2,0}; -- epsilon in the text
r={1,1,1,0}; -- theta in the text
s={0,2,1,1}; -- kappa in the text

S=ideal(x_(0,0,0)-t_1,x_(0,0,1)-t_2,x_(0,1,0)-t_3*c#0-t_5*d#0,x_(0,1,1),x_(0,2,0)-t_4*r#0-t_6*s#0,x_(0,2,1)-t_3*c#3-t_5*d#3,x_(1,0,0)-t_3*c#1-t_5*d#1,x_(1,0,1),x_(1,1,0)-t_1,x_(1,1,1)-t_2,x_(1,2,0),x_(1,2,1)-t_4*r#3-t_6*s#3,x_(2,0,0)-t_4*r#1-t_6*s#1,x_(2,0,1)-t_3*c#2-t_5*d#2,x_(2,1,0),x_(2,1,1)-t_4*r#2-t_6*s#2,x_(2,2,0)-t_1,x_(2,2,1)-t_2);

-- S is a parameterization of the 6-dimensional subspace (including coefficients of each basis element, t_1..t_6)

S'=gfanBuchberger(S);					   -- S' finds the reduced  Grobner basis of S (which is just row reduction in this case)

s'=drop(toList S' #1,-6);				   -- s' is S' with the terms containing t_1..t_6 removed (this means our desired 6-dimensional subspace is V(s'))

R=QQ[x_(0,0,0)..x_(2,2,1)];				   -- Redefine R to remove auxilliary variables t_1..t_6, and redefine S to be <s'>, an ideal of our newly defined ring R.
S=sub(ideal(s'),R);


col1=apply(flatten(table(3,2,(i,j)->x_(i,0,j)..x_(i,2,j))),toList); -- columns of 3x6 matrix formed by elements of the subspace
col2=apply(flatten(table(3,2,(i,j)->x_(0,i,j)..x_(2,i,j))),toList); -- columns of the partial transpose		

n1=transpose matrix col1;        
n2=transpose matrix col2;

D=minors(3,n1)+minors(3,n2);       -- Minors of 3x6 matrix and its partial transpose. D cuts out the tensors of border rank at most 2.
radical(S+D)			  
-- radical(S+D) cuts out the zero variety. Since the rationals are a perfect field, this implies that radical(S+D) also cuts out the zero variety
-- over the complexes. This exactly says that V(S) intersect V(D)={0}, i.e., S contains no non-zero tensors of border rank at most two.




-- We next verify that the 6-dimensional subspace of C^2 otimes C^2 otimes C^2 otimes C^2 presented in Equation (13) is 2-entangled. The methodology and
-- reasoning is exactly the same as in the previous case, so we avoid lengthy explanations here.

T=QQ[t_1..t_6,x_(0,0,0,0)..x_(1,1,1,1)];

Q=ideal(x_(0,0,0,0)-t_6,x_(0,0,0,1)-t_1,x_(0,1,0,0)-t_3-t_4,x_(0,1,0,1),x_(0,0,1,0)-t_2-t_4,x_(0,0,1,1)-t_5-2*t_6,x_(0,1,1,0)-t_1,x_(0,1,1,1)-t_4,x_(1,0,0,0)-t_3,x_(1,0,0,1)-t_2,x_(1,1,0,0)-t_5-t_6,x_(1,1,0,1)-t_1-t_3,x_(1,0,1,0),x_(1,0,1,1)-t_3-t_4,x_(1,1,1,0)-t_2,x_(1,1,1,1)-t_5);
Q'=gfanBuchberger(Q);
q'=drop(toList Q' #1,-6);
T=QQ[x_(0,0,0,0)..x_(1,1,1,1)];
Q=sub(ideal(q'),T);

dol1=apply(flatten(table(2,2,(i,j)->(x_(0,i,0,j)..x_(1,i,1,j)))),toList);       
dol2=apply(flatten(table(2,2,(i,j)->(x_(0,i,j,0)..x_(1,i,j,1)))),toList);
dol3=apply(flatten(table(2,2,(i,j)->(x_(i,j,0,0)..x_(i,j,1,1)))),toList);
-- These are the  columns of 4x4  matrix formed by elements of the subspace, and the other two flattenings.

m1=transpose matrix dol1;
m2=transpose matrix dol2;
m3=transpose matrix dol3;

J=minors(3,m1)+minors(3,m2)+minors(3,m3); -- J cuts out the variety of tensors of border rank at most 2.


radical(Q+J)          --radical(Q+J) cuts out the zero variety, so Q is 2-entangled.
