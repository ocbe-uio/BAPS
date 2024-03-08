function [muutokset, diffInCounts] = laskeMuutokset2( ...
  i1, globalRows, data, adjprior, priorTerm);
  % Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mikה olisi
  % muutos logml:ssה, mikהli korin i1 kaikki yksilצt siirretההn
  % koriin i.

  global COUNTS;      global SUMCOUNTS;
  global PARTITION;   global POP_LOGML;
  npops = size(COUNTS,3);
  muutokset = zeros(npops,1);

  i1_logml = POP_LOGML(i1);

  inds = find(PARTITION==i1);
  ninds = length(inds);

  if ninds==0
    diffInCounts = zeros(size(COUNTS,1), size(COUNTS,2));
    return;
  end

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
  i2_logml = POP_LOGML(i2);

  COUNTS(:,:,i2) = COUNTS(:,:,i2)+repmat(diffInCounts, [1 1 npops-1]);
  SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)+repmat(diffInSumCounts,[npops-1 1]);
  new_i2_logml = computePopulationLogml(i2, adjprior, priorTerm);
  COUNTS(:,:,i2) = COUNTS(:,:,i2)-repmat(diffInCounts, [1 1 npops-1]);
  SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)-repmat(diffInSumCounts,[npops-1 1]);

  muutokset(i2) = new_i1_logml - i1_logml ...
  + new_i2_logml - i2_logml;
end
