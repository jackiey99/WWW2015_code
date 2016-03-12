function score = label_fast_exact(aa,LL,currentTeam,i0, prune)
%proposed TEAMREP-FAST-EXACT for team replacement on labeled graph
%input:
%aa: the whole social network, e.g., author-author network
%L: label matrice cell, e.g., if there are dn skills, then L is a cell of size dn, 
%L{i} is a nxn diagonal matrix, L{i}(j,j) shows the strength of j-th person
%having i-th skill

%currentTeam: current members of a team, e.g., authors for a paper
%i0: the member to be replaced
%prune: prune or not?
%
%output:
%score: each row is a score and its candidate id, note it's not sorted

%Author: Liangyue Li (Arizona State University)

if nargin < 4
    prune = false;
end

n=size(aa,1);
dn=length(LL);

remainTeam = setdiff(currentTeam,i0);
currentTeam = [remainTeam, i0];

W = aa(currentTeam,currentTeam);
W = (triu(W,1) + tril(W,-1));

%W0 is the common part
n0 = length(currentTeam);
W0 = W;
W0(n0,1:n0)=0;
W0(1:n0,n0)=0;

L=cell(1,dn);
for i=1:dn
    L{i}=LL{i}(currentTeam,currentTeam);
end

%L0 is the common part
L0 = cell(1,dn);
for i=1:dn
    temp = L{i};
    temp(n0,n0)=0;
    L0{i}=temp;
end

cand = setdiff((1:n),currentTeam);

%Prune unpromising candidates
if prune == true
    cand = cand(sum(aa(cand,remainTeam),2)>0);
end

c=0.00000001;

q = ones(n0,1)/n0;
p = ones(n0,1)/n0;
qx = kron(q,q);
px = kron(p,p);

temp=zeros(n0*n0);
for i=1:dn
    temp = temp + kron(L{i}*W,L0{i}*W0);
end
invZ = inv(eye(n0*n0)-c*temp);

R =zeros(n0*n0,1);
for i=1:dn
    R = R + kron(L{i}*p,L0{i}*p);
end

base=qx'*invZ*R;
l = c*qx'*invZ;
r=invZ*R;

score = zeros(length(cand),2);

for i=1:length(cand)

    s=[zeros(1,n0-1),1];
    t=[aa(cand(i),remainTeam),0];
    A=[t',s'];
    B = [s;t];   
    
    a =[zeros(n0-1,1);1];
    b=cell(1,dn);
    for j=1:dn
        b{j}=[zeros(1,n0-1),LL{j}(cand(i),cand(i))];
    end
    
    A1=[];
    B1=[];
    for j=1:dn
        A1=[A1,kron(L{j},a)];
        B1=[B1;kron(eye(n0),b{j})]; 
    end
    
    X1=zeros(n0*n0,n0*2);
    X2=zeros(n0*n0,n0*2);
    for j=1:dn
        X1=X1+kron(L{j}*W,L0{j}*A);
        X2=X2+kron(L{j}*W,a*b{j}*A);
    end
    
    Y1=B1*kron(W,W0);
    Y2=kron(eye(n0),B);
    X=[A1,X1,X2];
    Y=[Y1;Y2;Y2];
    
    M=inv(eye((dn+4)*n0)-c*Y*invZ*X);
    r0=invZ*A1*B1*px;
  
    score(i,1) = base+qx'*r0+l*X*M*Y*(r+r0);
    score(i,2) = cand(i);
end


end
