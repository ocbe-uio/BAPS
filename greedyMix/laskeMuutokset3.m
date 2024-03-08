function muutokset = laskeMuutokset3(T2, inds2, globalRows, ...
  data, adjprior, priorTerm, i1)
  % Palauttaa length(unique(T2))*npops taulun, jossa (i,j):s alkio
  % kertoo, mikה olisi muutos logml:ssה, jos populaation i1 osapopulaatio
  % inds2(find(T2==i)) siirretההn koriin j.

  global COUNTS;
  global SUMCOUNTS;
  global PARTITION;
  global POP_LOGML;
  npops = size(COUNTS,3);
  npops2 = length(unique(T2));
  muutokset = zeros(npops2, npops);

  i1_logml = POP_LOGML(i1);
  for pop2 = 1:npops2
    inds = inds2(find(T2==pop2));
    ninds = length(inds);
    if ninds>0
      rows = [];
      for i = 1:ninds
        ind = inds(i);
        lisa = globalRows(ind,1):globalRows(ind,2);
        rows = [rows; lisa'];
        %rows = [rows; globalRows{ind}'];
      end
      diffInCounts = computeDiffInCounts(rows', size(COUNTS,1), size(COUNTS,2), data);
      diffInSumCounts = sum(diffInCounts);

      COUNTS(:,:,i1) = COUNTS(:,:,i1)-diffInCounts;
      SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)-diffInSumCounts;
      new_i1_logml = computePopulationLogml(i1, adjprior, priorTerm);
      COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
      SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;

      i2 = [1:i1-1 , i1+1:npops];
      i2_logml = POP_LOGML(i2)';

      COUNTS(:,:,i2) = COUNTS(:,:,i2)+repmat(diffInCounts, [1 1 npops-1]);
      SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)+repmat(diffInSumCounts,[npops-1 1]);
      new_i2_logml = computePopulationLogml(i2, adjprior, priorTerm)';
      COUNTS(:,:,i2) = COUNTS(:,:,i2)-repmat(diffInCounts, [1 1 npops-1]);
      SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)-repmat(diffInSumCounts,[npops-1 1]);

      muutokset(pop2,i2) = new_i1_logml - i1_logml ...
      + new_i2_logml - i2_logml;
    end
  end
end
