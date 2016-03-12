function ids = topfive(score)
    score = sortrows(score,-1);
    ids = score(1:5,2);
end