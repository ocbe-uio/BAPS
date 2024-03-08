function [muutokset, diffInCounts] = ...
  laskeMuutokset(ind, globalRows, data, adjprior, priorTerm)
  % Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mikה olisi
  % muutos logml:ssה, mikהli yksilצ ind siirretההn koriin i.
  % diffInCounts on poistettava COUNTS:in siivusta i1 ja lisהttהvה
  % COUNTS:in siivuun i2, mikהli muutos toteutetaan.
  %
  % Lisהys 25.9.2007:
  % Otettu kהyttצצn globaali muuttuja LOGDIFF, johon on tallennettu muutokset
  % logml:ssה siirrettהessה yksilצitה toisiin populaatioihin.

  global COUNTS;      global SUMCOUNTS;
  global PARTITION;   global POP_LOGML;
  global LOGDIFF;

  npops = size(COUNTS,3);
  muutokset = LOGDIFF(ind,:);

  i1 = PARTITION(ind);
  i1_logml = POP_LOGML(i1);
  muutokset(i1) = 0;

  rows = globalRows(ind,1):globalRows(ind,2);
  diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
  diffInSumCounts = sum(diffInCounts);

  COUNTS(:,:,i1) = COUNTS(:,:,i1)-diffInCounts;
  SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)-diffInSumCounts;
  new_i1_logml = computePopulationLogml(i1, adjprior, priorTerm);
  COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
  SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;

  i2 = find(muutokset==-Inf);     % Etsitההn populaatiot jotka muuttuneet viime kerran jהlkeen.
  i2 = setdiff(i2,i1);
  i2_logml = POP_LOGML(i2);

  ni2 = length(i2);

  COUNTS(:,:,i2) = COUNTS(:,:,i2)+repmat(diffInCounts, [1 1 ni2]);
  SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)+repmat(diffInSumCounts,[ni2 1]);
  new_i2_logml = computePopulationLogml(i2, adjprior, priorTerm);
  COUNTS(:,:,i2) = COUNTS(:,:,i2)-repmat(diffInCounts, [1 1 ni2]);
  SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)-repmat(diffInSumCounts,[ni2 1]);

  muutokset(i2) = new_i1_logml - i1_logml ...
  + new_i2_logml - i2_logml;
  LOGDIFF(ind,:) = muutokset;
end
