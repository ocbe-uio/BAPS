function popLogml = computePopulationLogml(pops, adjprior, priorTerm)
  % Palauttaa length(pops)*1 taulukon, jossa on laskettu korikohtaiset
  % logml:t koreille, jotka on mההritelty pops-muuttujalla.

  global COUNTS;
  global SUMCOUNTS;
  x = size(COUNTS,1);
  y = size(COUNTS,2);
  z = length(pops);

  popLogml = ...
  squeeze(sum(sum(reshape(...
  gammaln(repmat(adjprior,[1 1 length(pops)]) + COUNTS(:,:,pops)) ...
  ,[x y z]),1),2)) - sum(gammaln(1+SUMCOUNTS(pops,:)),2) - priorTerm;
end
