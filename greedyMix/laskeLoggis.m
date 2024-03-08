function loggis = laskeLoggis(counts, sumcounts, adjprior)
  npops = size(counts,3);

  logml2 = sum(sum(sum(gammaln(counts+repmat(adjprior,[1 1 npops]))))) ...
  - npops*sum(sum(gammaln(adjprior))) - ...
  sum(sum(gammaln(1+sumcounts)));
  loggis = logml2;
end
