%{
This script showcases how to use the proposed team member replacement
algorithms proposed in the WWW'15 paper:
Liangyue Li, Hanghang Tong, Nan Cao, Kate Ehrlich, Yu-Ru Lin, Norbou Buchler. Replacing the Irreplaceable: 
Fast Algorithm for Team Member Recommendation. International World Wide Web Conference (WWW), 2015.

Author: Liangyue Li (Arizona State University)

%}

load DBLP;

currentTeam = [916232, 250177, 219532, 545532, 756929];
i0= 250177;


dn=43; % number of skills
L = cell(1,dn);
for i=1:dn 
    L{i} = diag(count_label(:,i));
end

fileID=fopen('authorDict.txt');
authorDict=textscan(fileID,'%s','delimiter','\n');
authorDict=authorDict{1};

fprintf('We need to replace %s ...\n', authorDict{i0});

score = label_direct_recommend(aa,L,currentTeam,i0,true);
top5 = topfive(score);
display 'Using TEAMREP-BASIC after pruning, the top five candidates are:';
fprintf('%s \n', authorDict{top5});


score = label_fast_exact(aa,L,currentTeam,i0,true);
top5 = topfive(score);
display 'Using TEAMREP-FAST-EXACT, the top five candidates are:';
fprintf('%s \n', authorDict{top5});


score = label_fast_approx(aa,L,currentTeam,i0,true,2);
top5 = topfive(score);
display 'Using TEAMREP-FAST-APPROX, the top five candidates are:';
fprintf('%s \n', authorDict{top5});