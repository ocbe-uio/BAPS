function [newData, rowsFromInd, alleleCodes, noalle, adjprior, priorTerm] = ...
  handleData(raw_data)
  % Alkuperהisen datan viimeinen sarake kertoo, milt?yksilצlt?
  % kyseinen rivi on perהisin. Funktio tutkii ensin, ett?montako
  % rivi?maksimissaan on perהisin yhdelt?yksilצlt? jolloin saadaan
  % tietהה onko kyseess?haploidi, diploidi jne... Tהmהn jהlkeen funktio
  % lisהה tyhji?rivej?niille yksilצille, joilta on perהisin vהhemmהn
  % rivej?kuin maksimimההr?
  %   Mikהli jonkin alleelin koodi on =0, funktio muuttaa tהmהn alleelin
  % koodi pienimmהksi koodiksi, joka isompi kuin mikההn kהytצss?oleva koodi.
  % Tהmהn jהlkeen funktio muuttaa alleelikoodit siten, ett?yhden lokuksen j
  % koodit saavat arvoja vהlill?1,...,noalle(j).
  data = raw_data;
  nloci=size(raw_data,2)-1;

  dataApu = data(:,1:nloci);
  nollat = find(dataApu==0);
  if ~isempty(nollat)
    isoinAlleeli = max(max(dataApu));
    dataApu(nollat) = isoinAlleeli+1;
    data(:,1:nloci) = dataApu;
  end
  dataApu = []; nollat = []; isoinAlleeli = [];

  noalle=zeros(1,nloci);
  alleelitLokuksessa = cell(nloci,1);
  for i=1:nloci
    alleelitLokuksessaI = unique(data(:,i));
    alleelitLokuksessa{i,1} = alleelitLokuksessaI(find(alleelitLokuksessaI>=0));
    noalle(i) = length(alleelitLokuksessa{i,1});
  end
  alleleCodes = zeros(max(noalle),nloci);
  for i=1:nloci
    alleelitLokuksessaI = alleelitLokuksessa{i,1};
    puuttuvia = max(noalle)-length(alleelitLokuksessaI);
    alleleCodes(:,i) = [alleelitLokuksessaI; zeros(puuttuvia,1)];
  end

  for loc = 1:nloci
    for all = 1:noalle(loc)
      data(find(data(:,loc)==alleleCodes(all,loc)), loc)=all;
    end;
  end;

  nind = max(data(:,end));
  nrows = size(data,1);
  ncols = size(data,2);
  rowsFromInd = zeros(nind,1);
  for i=1:nind
    rowsFromInd(i) = length(find(data(:,end)==i));
  end
  maxRowsFromInd = max(rowsFromInd);
  a = -999;
  emptyRow = repmat(a, 1, ncols);
  lessThanMax = find(rowsFromInd < maxRowsFromInd);
  missingRows = maxRowsFromInd*nind - nrows;
  data = [data; zeros(missingRows, ncols)];
  pointer = 1;
  for ind=lessThanMax'    %Kהy lהpi ne yksilצt, joilta puuttuu rivej?
    miss = maxRowsFromInd-rowsFromInd(ind);  % Tהlt?yksilצlt?puuttuvien lkm.
    for j=1:miss
      rowToBeAdded = emptyRow;
      rowToBeAdded(end) = ind;
      data(nrows+pointer, :) = rowToBeAdded;
      pointer = pointer+1;
    end
  end
  data = sortrows(data, ncols);   % Sorttaa yksilצiden mukaisesti
  newData = data;
  rowsFromInd = maxRowsFromInd;

  adjprior = zeros(max(noalle),nloci);
  priorTerm = 0;
  for j=1:nloci
    adjprior(:,j) = [repmat(1/noalle(j), [noalle(j),1]) ; ones(max(noalle)-noalle(j),1)];
    priorTerm = priorTerm + noalle(j)*gammaln(1/noalle(j));
  end
end
