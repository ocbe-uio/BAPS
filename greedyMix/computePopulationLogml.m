function popLogml = computePopulationLogml(pops, adjprior, priorTerm)
  % Palauttaa length(pops)*1 taulukon, jossa on laskettu korikohtaiset
  % logml:t koreille, jotka on mההritelty pops-muuttujalla.

  global COUNTS;
  global SUMCOUNTS;
  x = size(COUNTS, 1);
  y = size(COUNTS, 2);
  z = length(pops);

  rep_adj = repmat(adjprior, [1 1 z]);
  gamma_rep_counts = gammaln(rep_adj + COUNTS(:, :, pops));
  gamma_sum_counts = sum(gammaln(1 + SUMCOUNTS(pops, :)), 2);
  gamma_rep_counts_sum = sum(sum(reshape(gamma_rep_counts, [x y z]), 1), 2);
  gamma_rep_counts_reshaped = squeeze(gamma_rep_counts_sum);
  popLogml = gamma_rep_counts_reshaped - gamma_sum_counts - priorTerm;
end
