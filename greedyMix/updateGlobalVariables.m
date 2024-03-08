function updateGlobalVariables(ind, i2, diffInCounts, ...
  adjprior, priorTerm)
  % Suorittaa globaalien muuttujien muutokset, kun yksilצ ind
  % on siirretההn koriin i2.

  global PARTITION;
  global COUNTS;
  global SUMCOUNTS;
  global POP_LOGML;
  global LOGDIFF;

  i1 = PARTITION(ind);
  PARTITION(ind)=i2;

  COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
  COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
  SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - sum(diffInCounts);
  SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + sum(diffInCounts);

  POP_LOGML([i1 i2]) = computePopulationLogml([i1 i2], adjprior, priorTerm);

  LOGDIFF(:,[i1 i2]) = -Inf;
  inx = [find(PARTITION==i1); find(PARTITION==i2)];
  LOGDIFF(inx,:) = -Inf;
end
