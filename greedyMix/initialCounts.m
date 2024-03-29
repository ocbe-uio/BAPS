function [sumcounts, counts, logml] = ...
  initialCounts(partition, data, npops, rows, noalle, adjprior)

  nloci=size(data,2);
  ninds = size(rows, 1);

  koot = rows(:,1) - rows(:,2) + 1;
  maxSize = max(koot);

  counts = zeros(max(noalle),nloci,npops);
  sumcounts = zeros(npops,nloci);
  for i=1:npops
    for j=1:nloci
      havainnotLokuksessa = find(partition==i & data(:,j)>=0);
      sumcounts(i,j) = length(havainnotLokuksessa);
      for k=1:noalle(j)
        alleleCode = k;
        N_ijk = length(find(data(havainnotLokuksessa,j)==alleleCode));
        counts(k,j,i) = N_ijk;
      end
    end
  end

  %initializeGammaln(ninds, maxSize, max(noalle));

  logml = laskeLoggis(counts, sumcounts, adjprior);
end
