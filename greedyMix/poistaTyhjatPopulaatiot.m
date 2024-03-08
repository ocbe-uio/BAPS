function npops = poistaTyhjatPopulaatiot(npops)
  % Poistaa tyhjentyneet populaatiot COUNTS:ista ja
  % SUMCOUNTS:ista. Pהivittהה npops:in ja PARTITION:in.

  global COUNTS;
  global SUMCOUNTS;
  global PARTITION;
  global LOGDIFF;

  notEmpty = find(any(SUMCOUNTS,2));
  COUNTS = COUNTS(:,:,notEmpty);
  SUMCOUNTS = SUMCOUNTS(notEmpty,:);
  LOGDIFF = LOGDIFF(:,notEmpty);

  for n=1:length(notEmpty)
    apu = find(PARTITION==notEmpty(n));
    PARTITION(apu)=n;
  end
  npops = length(notEmpty);
end
