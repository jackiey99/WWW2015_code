function score = label_direct_recommend(aa,L,currentTeam,i0,prune)
%The proposed TEAMREP-BASIC algorithm
%input:
%aa: the whole social network, e.g., author-author network
%L: label matrice cell, e.g., if there are dn skills, then L is a cell of size dn, 
%L{i} is a nxn diagonal matrix, L{i}(j,j) shows the strength of j-th person
%having i-th skill
%
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
remainTeam = setdiff(currentTeam,i0);
currentTeam = [remainTeam, i0];

A1 = aa(currentTeam,currentTeam);
A1 = (triu(A1,1) + tril(A1,-1));  % remove diagonal elements

cand = setdiff((1:n),currentTeam); % this is the set of candidates

% prune those unpromising candidates
if prune == true
    cand = cand(sum(aa(cand,remainTeam),2)>0);
end

% decay factor in RWR, set it so that it converges
c=0.00000001;

% number of skills
dn=length(L);

score = zeros(length(cand), 2);

for i=1:length(cand)
    newTeam = [remainTeam,cand(i)];
    LL = zeros(length(newTeam)^2);
    for j=1:dn
        LL = LL + kron(L{j}(currentTeam,currentTeam),L{j}(newTeam,newTeam));
    end
    A2 = aa(newTeam,newTeam);
    A2 = (triu(A2,1) + tril(A2,-1));
   
    score(i,1) = label_gs(A1,A2,LL,c);
    score(i,2) = cand(i);
end

end

function sim = label_gs(A,B,L,c)
% graph kernel computation for attributes graphs
% input:
%   A: first graph
%   B: second graph
%   L: label matrix
%   c: decay factor
% output:
%   sim: similarity between A and B

    n1 = size(A,1);
    n2 = size(B,1);
    q = {ones(n1,1)/n1,ones(n2,1)/n2};
    p = {ones(n1,1)/n1,ones(n2,1)/n2};
    p1 = p{1};
    p2 = p{2};
    q1 = q{1};
    q2 = q{2};
    
    X = kron(A,B);
    qx = kron(q1,q2);
    px = kron(p1,p2);
    sim = qx' * inv(eye(n1 * n2) - c * L * X) * L * px;
end