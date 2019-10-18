function [ Kb, lambda, v11, v11inv, Mb, p, strnMode] = getSBFEMStiffMat_2NodeEle( xy, eConn, D, density, g)
% solution of a polygon element by SBFEM

% xy(i,:)       x-, y-coordinates of node i
% eConn(i,:)    nodes of 2-node element i
% D             elasticity matrix 3x3
% density 
% g             gravitational acceleration

%% coefficient matrices
nn = size(eConn,2); % number of line element
nd = nn+nn;         % number of degrees of freedom (DOF)
id = 1:nd;          % working vector
E0 = zeros(nd, nd); E1 = E0; E2 = E0; M0 = E0; % initialisation
n1 = eConn(1,:); %1st node of elements
n2 = eConn(2,:); %2nd node of elements 
dof = [ n1+n1-1; n1+n1; n2+n2-1; n2+n2 ]; %DOFs 
x = xy(:,1)'; dx = (x(n2)-x(n1))/2; ax = (x(n2)+x(n1))/2;
y = xy(:,2)'; dy = (y(n2)-y(n1))/2; ay = (y(n2)+y(n1))/2;
J2 = 2*(ax.*dy-ay.*dx);
for ie = 1:nn
    C1 = [ dy(ie) 0;   0 -dx(ie);   -dx(ie) dy(ie) ];
    C2 = [ ay(ie) 0;   0 -ax(ie);   -ax(ie) ay(ie) ];
    a = 1/J2(ie); Q0 = a/3*(C1'*D*C1); Q1 = -a*(C2'*D*C1); Q2 = a*(C2'*D*C2);
    d = dof(:,ie);
    E0(d,d)=E0(d,d)+[ 4*Q0 2*Q0; 2*Q0 4*Q0];
    E1(d,d)=E1(d,d)+[ -Q0-Q1  Q0-Q1;  Q0+Q1 -Q0+Q1 ];
    E2(d,d)=E2(d,d)+[  Q0+Q2 -Q0-Q2; -Q0-Q2  Q0+Q2 ];
    b = J2(ie)*density/6;
    M0(d,d)=M0(d,d)+[ b+b 0  b 0; 0 b+b  0 b; b 0  b+b 0; 0 b 0 b+b ];
end
%% stiffness matrix
m = E0\[E1' -eye(nd)]; Z = [m; E1*m(:,id)-E2 -m(:,id)'];
[v, d] = eig(Z);  lambda = diag(d);
[~, idx] = sort(real(lambda),'ascend'); %sort eignvalues in ascending order
lambda = lambda(idx(id))';  v = v(:, idx(id)); %rearrange eigenvalues and eigenvectors
lambda(end-1:end) = 0;  v(:,end-1:end) = 0;%rigid body translational modes
v(1:2:nd,end-1) = 1;  v(2:2:nd,end) = 1;
v11 = v(id, :); v11inv = inv(v11);
Kb  = real(v(nd+id, :)*v11inv);%stiffness matrix

%% mass matrix; 
M0 = v11'*M0*v11;
am = lambda(ones(1,nd),:); 
M0 = M0./(2-am-am'); 
Mb = real(v11inv'*M0*v11inv);
%% self-weight
p = real(g*Mb*v11(:,end));

%% strain modes
v11b = -v11(:,1:end-2).*am(:,1:end-2);
strnMode = zeros(3*nn, nd-2);
for ie = 1:nn
    au = (v11b(dof(3:4,ie),:)+v11b(dof(1:2,ie),:));
    du = (v11(dof(3:4,ie),1:end-2)-v11(dof(1:2,ie),1:end-2));
    b1h = [ dy(ie) 0;   0 -dx(ie);   -dx(ie) dy(ie) ]/J2(ie);
    b2h = [ -ay(ie) 0;   0 ax(ie);   ax(ie) -ay(ie) ]/J2(ie);
    strnMode(3*(ie-1)+1:3*ie,:) = b1h*au + b2h*du;
end

