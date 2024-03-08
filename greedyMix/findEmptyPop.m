function [emptyPop, pops] = findEmptyPop(npops)
  % Palauttaa ensimmהisen tyhjהn populaation indeksin. Jos tyhjiה
  % populaatioita ei ole, palauttaa -1:n.

  global PARTITION;
  pops = unique(PARTITION)';
  if (length(pops) ==npops)
    emptyPop = -1;
  else
    popDiff = diff([0 pops npops+1]);
    emptyPop = min(find(popDiff > 1));
  end
end
