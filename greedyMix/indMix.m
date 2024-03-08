function [logml, npops, partitionSummary] = indMix(c, npops, dispText)
  % Greedy search algorithm with unknown number of classes for regular
  % clustering.
  % Input npops is not used if called by greedyMix or greedyPopMix.

  logml = 1;

  global PARTITION;
  global COUNTS;
  global SUMCOUNTS;
  global POP_LOGML;
  global LOGDIFF;
  clearGlobalVars;

  noalle = c.noalle;
  rows = c.rows;
  data = c.data;
  adjprior = c.adjprior;
  priorTerm = c.priorTerm;
  rowsFromInd = c.rowsFromInd;

  if isfield(c,'dist')
    dist = c.dist; Z = c.Z;
  end

  clear c;

  if nargin < 2
    dispText = 1;
    npopstext = [];
    ready = false;
    teksti = 'Input upper bound to the number of populations (possibly multiple values): ';
    while ready == false
      npopstextExtra = inputdlg(teksti ,...
      'Input maximum number of populations',1,{'3'});
      drawnow
      if isempty(npopstextExtra)  % Painettu Cancel:ia
        return
      end
      npopstextExtra = npopstextExtra{1};
      if length(npopstextExtra)>=255
        npopstextExtra = npopstextExtra(1:255);
        npopstext = [npopstext ' ' npopstextExtra];
        teksti = 'The input field length limit (255 characters) was reached. Input more values: ';
      else
        npopstext = [npopstext ' ' npopstextExtra];
        ready = true;
      end
    end
    clear ready; clear teksti;
    if isempty(npopstext) | length(npopstext)==1
      return
    else
      npopsTaulu = str2num(npopstext);
      ykkoset = find(npopsTaulu==1);
      npopsTaulu(ykkoset) = [];   % Mikäli ykkösiä annettu ylärajaksi, ne poistetaan.
      if isempty(npopsTaulu)
        logml = 1; partitionSummary=1; npops=1;
        return
      end
      clear ykkoset;
    end
  else
    npopsTaulu = npops;
  end

  nruns = length(npopsTaulu);

  initData = data;
  data = data(:,1:end-1);

  logmlBest = -1e50;
  partitionSummary = -1e50*ones(30,2);  % Tiedot 30 parhaasta partitiosta (npops ja logml)
  partitionSummary(:,1) = zeros(30,1);
  worstLogml = -1e50; worstIndex = 1;

  for run = 1:nruns
    npops = npopsTaulu(run);
    if dispText
      dispLine;
      disp(['Run ' num2str(run) '/' num2str(nruns) ...
      ', maximum number of populations ' num2str(npops) '.']);
    end
    ninds = size(rows,1);

    initialPartition = admixture_initialization(initData, npops, Z);
    [sumcounts, counts, logml] = ...
    initialCounts(initialPartition, data, npops, rows, noalle, adjprior);
    PARTITION = zeros(ninds, 1);
    for i=1:ninds
      apu = rows(i);
      PARTITION(i) = initialPartition(apu(1));
    end

    COUNTS = counts;
    SUMCOUNTS = sumcounts;
    POP_LOGML = computePopulationLogml(1:npops, adjprior, priorTerm);
    LOGDIFF = repmat(-Inf,ninds,npops);
    clear initialPartition; clear counts; clear sumcounts;

    % PARHAAN MIXTURE-PARTITION ETSIMINEN
    nRoundTypes = 7;
    kokeiltu = zeros(nRoundTypes, 1);
    roundTypes = [1 1];  %Ykkösvaiheen sykli kahteen kertaan.
    ready = 0; vaihe = 1;

    if dispText
      disp(' ');
      disp(['Mixture analysis started with initial ' num2str(npops) ' populations.']);
    end

    while ready ~= 1
      muutoksia = 0;

      if dispText
        disp(['Performing steps: ' num2str(roundTypes)]);
      end

      for n = 1:length(roundTypes)

        round = roundTypes(n);
        kivaluku=0;

        if kokeiltu(round) == 1     %Askelta kokeiltu viime muutoksen jälkeen
          % do nothing
        elseif round==0 | round==1   %Yksilön siirtäminen toiseen populaatioon.
          inds = 1:ninds;
          aputaulu = [inds' rand(ninds,1)];
          aputaulu = sortrows(aputaulu,2);
          inds = aputaulu(:,1)';
          muutosNyt = 0;
          for ind = inds
            i1 = PARTITION(ind);
            [muutokset, diffInCounts] = laskeMuutokset(ind, rows, ...
            data, adjprior, priorTerm);

            if round==1
              [maxMuutos, i2] = max(muutokset);
            end

            if (i1~=i2 & maxMuutos>1e-5)
              % Tapahtui muutos
              muutoksia = 1;
              if muutosNyt == 0
                muutosNyt = 1;
                if dispText
                  disp('Action 1');
                end
              end
              kokeiltu = zeros(nRoundTypes,1);
              kivaluku = kivaluku+1;
              updateGlobalVariables(ind, i2, diffInCounts,...
              adjprior, priorTerm);
              logml = logml+maxMuutos;
              if logml>worstLogml
                [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex);
                if (added==1)
                  [worstLogml, worstIndex] = min(partitionSummary(:,2));
                end
              end
            end
          end
          if muutosNyt == 0
            kokeiltu(round) = 1;
          end

        elseif round==2  %Populaation yhdistäminen toiseen.
          maxMuutos = 0;
          for pop = 1:npops
            [muutokset, diffInCounts] = laskeMuutokset2(pop, rows, ...
            data, adjprior, priorTerm);
            [isoin, indeksi] = max(muutokset);
            if isoin>maxMuutos
              maxMuutos = isoin;
              i1 = pop;
              i2 = indeksi;
              diffInCountsBest = diffInCounts;
            end
          end

          if maxMuutos>1e-5
            muutoksia = 1;
            kokeiltu = zeros(nRoundTypes,1);
            updateGlobalVariables2(i1,i2, diffInCountsBest, ...
            adjprior, priorTerm);
            logml = logml + maxMuutos;
            if dispText
              disp('Action 2');
            end
            if logml>worstLogml
              [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex);
              if (added==1)
                [worstLogml, worstIndex] = min(partitionSummary(:,2));
              end
            end
          else
            kokeiltu(round) = 1;
          end


        elseif round==3 || round==4 %Populaation jakaminen osiin.
          maxMuutos = 0;
          ninds = size(rows,1);
          for pop = 1:npops
            inds2 = find(PARTITION==pop);
            ninds2 = length(inds2);
            if ninds2>2
              dist2 = laskeOsaDist(inds2, dist, ninds);
              Z2 = linkage(dist2');
              if round==3
                npops2 = max(min(20, floor(ninds2/5)),2);
              elseif round==4
                npops2 = 2;  %Moneenko osaan jaetaan
              end
              T2 = cluster_own(Z2, npops2);
              muutokset = laskeMuutokset3(T2, inds2, rows, data, ...
              adjprior, priorTerm, pop);
              [isoin, indeksi] = max(muutokset(1:end));
              if isoin>maxMuutos
                maxMuutos = isoin;
                muuttuvaPop2 = rem(indeksi,npops2);
                if muuttuvaPop2==0
                  muuttuvaPop2 = npops2;
                end
                muuttuvat = inds2(find(T2==muuttuvaPop2));
                i2 = ceil(indeksi/npops2);
              end
            end
          end
          if maxMuutos>1e-5
            muutoksia = 1;
            kokeiltu = zeros(nRoundTypes,1);
            rivit = [];
            for i = 1:length(muuttuvat)
              ind = muuttuvat(i);
              lisa = rows(ind,1):rows(ind,2);
              rivit = [rivit; lisa'];
            end
            diffInCounts = computeDiffInCounts(rivit', size(COUNTS,1), ...
            size(COUNTS,2), data);
            i1 = PARTITION(muuttuvat(1));
            updateGlobalVariables3(muuttuvat, diffInCounts, ...
            adjprior, priorTerm, i2);
            logml = logml + maxMuutos;
            if dispText
              if round==3
                disp('Action 3');
              else
                disp('Action 4');
              end
            end
            if logml>worstLogml
              [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex);
              if (added==1)
                [worstLogml, worstIndex] = min(partitionSummary(:,2));
              end
            end

          else
            kokeiltu(round)=1;
          end
        elseif round == 5 || round == 6
          j=0;
          muutettu = 0;
          poplogml = POP_LOGML;
          partition = PARTITION;
          counts = COUNTS;
          sumcounts = SUMCOUNTS;
          logdiff = LOGDIFF;

          pops = randperm(npops);
          while (j < npops & muutettu == 0)
            j = j+1;
            pop = pops(j);
            totalMuutos = 0;
            inds = find(PARTITION==pop);
            if round == 5
              aputaulu = [inds rand(length(inds),1)];
              aputaulu = sortrows(aputaulu,2);
              inds = aputaulu(:,1)';
            elseif round == 6
              inds = returnInOrder(inds, pop, rows, data, adjprior, priorTerm);
            end

            i = 0;

            while (length(inds) > 0 & i < length(inds))
              i = i+1;
              ind =inds(i);

              [muutokset, diffInCounts] = laskeMuutokset(ind, rows, ...
              data, adjprior, priorTerm);
              muutokset(pop) = -1e50;   % Varmasti ei suurin!!!
              [maxMuutos, i2] = max(muutokset);
              updateGlobalVariables(ind, i2, diffInCounts,...
              adjprior, priorTerm);

              totalMuutos = totalMuutos+maxMuutos;
              logml = logml+maxMuutos;
              if round == 6
                % Lopetetaan heti kun muutos on positiivinen.
                if totalMuutos > 1e-5
                  i=length(inds);
                end
              end
            end

            if totalMuutos>1e-5
              kokeiltu = zeros(nRoundTypes,1);
              muutettu=1;
              if muutoksia==0
                muutoksia = 1;  % Ulompi kirjanpito.
                if dispText
                  if round==5
                    disp('Action 5');
                  else
                    disp('Action 6');
                  end
                end
              end
              if logml>worstLogml
                [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex);
                if (added==1)
                  [worstLogml, worstIndex] = min(partitionSummary(:,2));
                end
              end
            else
              % Missään vaiheessa tila ei parantunut.
              % Perutaan kaikki muutokset.
              PARTITION = partition;
              SUMCOUNTS = sumcounts;
              POP_LOGML = poplogml;
              COUNTS = counts;
              logml = logml - totalMuutos;
              LOGDIFF = logdiff;
              kokeiltu(round)=1;
            end
          end
          clear partition; clear sumcounts; clear counts; clear poplogml;

        elseif round == 7
          emptyPop = findEmptyPop(npops);
          j = 0;
          pops = randperm(npops);
          muutoksiaNyt = 0;
          if emptyPop == -1
            j = npops;
          end
          while (j < npops)
            j = j +1;
            pop = pops(j);
            inds2 = find(PARTITION == pop);
            ninds2 = length(inds2);
            if ninds2 > 5
              partition = PARTITION;
              sumcounts = SUMCOUNTS;
              counts = COUNTS;
              poplogml = POP_LOGML;
              logdiff = LOGDIFF;

              dist2 = laskeOsaDist(inds2, dist, ninds);
              Z2 = linkage(dist2');
              T2 = cluster_own(Z2, 2);
              muuttuvat = inds2(find(T2 == 1));

              muutokset = laskeMuutokset3(T2, inds2, rows, data, ...
              adjprior, priorTerm, pop);
              totalMuutos = muutokset(1, emptyPop);

              rivit = [];
              for i = 1:length(muuttuvat)
                ind = muuttuvat(i);
                lisa = rows(ind,1):rows(ind,2);
                rivit = [rivit lisa];
              end
              diffInCounts = computeDiffInCounts(rivit, size(COUNTS,1), ...
              size(COUNTS,2), data);

              updateGlobalVariables3(muuttuvat, diffInCounts, ...
              adjprior, priorTerm, emptyPop);

              muutettu = 1;
              while (muutettu == 1)
                muutettu = 0;
                % Siirretään yksilöitä populaatioiden välillä
                muutokset = laskeMuutokset5(inds2, rows, data, ...
                adjprior, priorTerm, pop, emptyPop);

                [maxMuutos, indeksi] = max(muutokset);

                muuttuva = inds2(indeksi);
                if (PARTITION(muuttuva) == pop)
                  i2 = emptyPop;
                else
                  i2 = pop;
                end

                if maxMuutos > 1e-5
                  rivit = rows(muuttuva,1):rows(muuttuva,2);
                  diffInCounts = computeDiffInCounts(rivit, size(COUNTS,1), ...
                  size(COUNTS,2), data);
                  updateGlobalVariables3(muuttuva,diffInCounts, ...
                  adjprior, priorTerm, i2);
                  muutettu = 1;
                  totalMuutos = totalMuutos + maxMuutos;
                end

              end

              if totalMuutos > 1e-5
                muutoksia = 1;
                logml = logml + totalMuutos;
                if logml>worstLogml
                  [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex);
                  if (added==1)
                    [worstLogml, worstIndex] = min(partitionSummary(:,2));
                  end
                end
                if muutoksiaNyt == 0
                  if dispText
                    disp('Action 7');
                  end
                  muutoksiaNyt = 1;
                end
                kokeiltu = zeros(nRoundTypes,1);
                j = npops;
              else
                %palutetaan vanhat arvot
                PARTITION = partition;
                SUMCOUNTS = sumcounts;
                COUNTS = counts;
                POP_LOGML = poplogml;
                LOGDIFF = logdiff;
              end

            end

          end

          if muutoksiaNyt == 0
            kokeiltu(round)=1;
          end

        end
      end

      if muutoksia == 0
        if vaihe==1
          vaihe = 2;
        elseif vaihe==2
          vaihe = 3;
        elseif vaihe==3
          vaihe = 4;
        elseif vaihe==4;
          vaihe = 5;
        elseif vaihe==5
          ready = 1;
        end
      else
        muutoksia = 0;
      end

      if ready==0
        if vaihe==1
          roundTypes=[1];
        elseif vaihe==2
          roundTypes = [2 1];
        elseif vaihe==3
          roundTypes=[5 5 7];
        elseif vaihe==4
          roundTypes=[4 3 1];
        elseif vaihe==5
          roundTypes=[6 7 2 3 4 1];
        end
      end
    end

    % TALLENNETAAN

    npops = poistaTyhjatPopulaatiot(npops);
    POP_LOGML = computePopulationLogml(1:npops, adjprior, priorTerm);
    if dispText
      disp(['Found partition with ' num2str(npops) ' populations.']);
      disp(['Log(ml) = ' num2str(logml)]);
      disp(' ');
    end

    if logml>logmlBest
      % Päivitetään parasta löydettyä partitiota.
      logmlBest = logml;
      npopsBest = npops;
      partitionBest = PARTITION;
      countsBest = COUNTS;
      sumCountsBest = SUMCOUNTS;
      pop_logmlBest = POP_LOGML;
      logdiffbest = LOGDIFF;
    end

  end

  logml = logmlBest;
  npops = npopsBest;
  PARTITION = partitionBest;
  COUNTS = countsBest;
  SUMCOUNTS = sumCountsBest;
  POP_LOGML = pop_logmlBest;
  LOGDIFF = logdiffbest;
end
