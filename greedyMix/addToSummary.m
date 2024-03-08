function [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex)
  % Tiedetההn, ettה annettu logml on isompi kuin huonoin arvo
  % partitionSummary taulukossa. Jos partitionSummary:ssה ei vielה ole
  % annettua logml arvoa, niin lisהtההn worstIndex:in kohtaan uusi logml ja
  % nykyistה partitiota vastaava nclusters:in arvo. Muutoin ei tehdה mitההn.

  apu = find(abs(partitionSummary(:,2)-logml)<1e-5);
  if isempty(apu)
    % Nyt lצydetty partitio ei ole vielה kirjattuna summaryyn.
    global PARTITION;
    npops = length(unique(PARTITION));
    partitionSummary(worstIndex,1) = npops;
    partitionSummary(worstIndex,2) = logml;
    added = 1;
  else
    added = 0;
  end
end
