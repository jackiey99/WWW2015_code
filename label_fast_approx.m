function score = label_fast_approx(aa,LL,currentTeam,i0,prune,r)
%proposed TEAMREP-FAST-APPROX for team replacement on labeled graph
%input:
%aa: the whole social network, e.g., author-author network
%L: label matrice cell, e.g., if there are dn skills, then L is a cell of size dn, 
%L{i} is a nxn diagonal matrix, L{i}(j,j) shows the strength of j-th person
%having i-th skill
%
%currentTeam: current members of a team, e.g., authors for a paper
%i0: the member to be replaced
%prune: prune or not?
%r: the approximate rank 
%
%output:
%score: each row is a score and its candidate id, note it's not sorted

%Author: Liangyue Li (Arizona State University)


if nargin < 4
    prune = false;
end

n=size(aa,1);
remainTeam = setdiff(currentTeam,i0);
currentTeam = [remainTeam, i0];

n0 = length(currentTeam);
W0 = zeros(n0,n0);
w_temp = aa(remainTeam,remainTeam);
w_temp = (triu(w_temp,1) + tril(w_temp,-1));
W0(1:n0-1,1:n0-1) = w_temp;

if nargin < 6
    r = floor(sqrt(n0));
end

% top r eigen-decomposition
[U0,Lam0,~]=svds(W0,r,'L');
U = U0;
V = Lam0*U0';
s = [zeros(n0-1,1);1];

w_original = aa(i0,currentTeam);
X = [U,w_original', s];
Y = [V;s';w_original];

cand = setdiff((1:n),currentTeam);
if prune == true
    cand = cand(sum(aa(cand,remainTeam),2)>0);
end

c=0.00000001;

dn=length(LL);

q = {ones(n0,1)/n0,ones(n0,1)/n0};
p = {ones(n0,1)/n0,ones(n0,1)/n0};
p1 = p{1};
p2 = p{2};
q1 = q{1};
q2 = q{2};

score = zeros(length(cand), 2);

for i=1:length(cand)
    newTeam = [remainTeam,cand(i)];
    w_new = aa(cand(i),newTeam);
    X_new = [U, w_new', s];
    Y_new = [V;s';w_new];

    for j=1:dn
        L1{j} = LL{j}(currentTeam,currentTeam);
        L2{j} = LL{j}(newTeam,newTeam);
    end
 
    temp  = zeros((r+2)^2,(r+2)^2);
    for j=1:dn
        temp = temp + kron(Y*L1{j}*X,Y_new*L2{j}*X_new);
    end
    
    tildeLam = inv( eye((r+2)^2) - c*temp);

    
    L = zeros(1,(r+2)^2);
    R = zeros( (r+2)^2,1);
    base = 0;
    for j=1:dn
        L = L + kron(q1'*L1{j}*X,q2'*L2{j}*X_new);
        R = R + kron(Y*L1{j}*p1,Y_new*L2{j}*p2);
        base = base + (q1'*L1{j}*p1)*(q2'*L2{j}*p2);
    end

    score(i,1) = base + c*L*tildeLam*R;
    score(i,2) = cand(i);
end

end