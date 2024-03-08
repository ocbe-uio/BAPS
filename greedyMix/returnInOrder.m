function inds = returnInOrder(inds, pop, globalRows, data, ...
  adjprior, priorTerm)
  % Palauttaa yksilצt jהrjestyksessה siten, ettה ensimmהisenה on
  % se, jonka poistaminen populaatiosta pop nostaisi logml:n
  % arvoa eniten.

  global COUNTS;
  global SUMCOUNTS;
  ninds = length(inds);
  apuTaulu = [inds, zeros(ninds,1)];

  for i=1:ninds
    ind =inds(i);
    rows = globalRows(i,1):globalRows(i,2);
    diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
    diffInSumCounts = sum(diffInCounts);

    diffCounts = COUNTS(:,:,pop) - diffInCounts; % TODO: workaround. Check if appropriate!
    COUNTS(:,:,pop) = max(diffCounts, 0); % Ensure non-negative values % TODO: workaround. Check if appropriate!
    SUMCOUNTS(pop,:) = SUMCOUNTS(pop,:) - diffInSumCounts;
    apuTaulu(i, 2) = computePopulationLogml(pop, adjprior, priorTerm);
    COUNTS(:,:,pop) = COUNTS(:,:,pop) + diffInCounts;
    SUMCOUNTS(pop,:) = SUMCOUNTS(pop,:) + diffInSumCounts;
  end
  apuTaulu = sortrows(apuTaulu,2);
  inds = apuTaulu(ninds:-1:1,1);
end
