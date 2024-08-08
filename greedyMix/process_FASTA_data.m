function processed_data = process_FASTA_data(file, partitionCompare, coordinates)
  if ~isempty(partitionCompare)
    fprintf(1, 'Data: %s\n', file);
  end
  [heds, seqs] = fastaread(file);
  seqs = seqs(:);
  alnMat = cell2mat(seqs);
  nSeq = length(seqs);
  clear seqs;

  cc = preprocAln(alnMat);
  cc.heds = heds;
  cc.nSeq = nSeq;
  dist = seqpdist(alnMat,'method','p-distance');
  Z = linkage(dist,'complete');

  ninds = nSeq;

  clear alnMat heds nSeq

  filename2 = coordinates;
  if isempty(filename2)
    return
  end

  coordinates = load(coordinates);
  [viallinen coordinates] = testaaKoordinaatit(ninds, coordinates); % added by Lu Cheng, 05.12.2012
  if viallinen
    disp('Incorrect coordinates');
    return
  end

  inp = [file ' & ' filename2];
  h0 = findobj('Tag','file_text');
  set(h0,'String',inp);
  clear h0; clear inp;
  clear file; clear filename2; clear pathname1; clear pathname2;

  input_pops = input(['When using data which are in FASTA-format, '...
  'you can specify the sampling populations of the individuals by '...
  'giving two additional files: one containing the names of the '...
  'populations, the other containing the indices of the first '...
  'individuals of the populations. Do you wish to specify the '...
  'sampling populations? (y/N)'], 's');
  if isequal(input_pops,'y')
    [namefile, namepath] = uigetfile('*.txt', 'Load population names');
    if namefile==0
      kysyToinen = 0;
    else
      kysyToinen = 1;
    end
    if kysyToinen==1
      [indicesfile, indicespath] = uigetfile('*.txt', 'Load population indices');
      if indicesfile==0
        popnames = [];
      else
        popnames = initPopNames([namepath namefile],[indicespath indicesfile]);
      end
    else
      popnames = [];
    end
  else
    popnames = [];
  end

  disp('Pre-processing the data. This may take several minutes.');

  [cliques, separators, vorPoints, vorCells, pointers] = ...
  handleCoords(coordinates);

  cc.locCliques = cliques;
  cc.locSeparators = separators;
  cc.popnames = popnames;
  cc.vorPoints = vorPoints;
  cc.vorCells = vorCells;
  cc.pointers = pointers;
  cc.coordinates = coordinates;
  format_type = 'FASTA';

  save_preproc = input('Do you wish to save pre-processed data? (y/N)', 's');
  if isequal(save_preproc,'y')
    [filename, pathname] = uiputfile('*.mat','Save pre-processed data as');
    kokonimi = [pathname filename];
    save(kokonimi,'cc','dist','Z','format_type','-v7.3'); % added by Lu Cheng, 08.06.2012
  end

  handleIndiFastaCase(cc,dist,Z);

end

function [viallinen coordinates] = testaaKoordinaatit(ninds, coordinates)
  % Testaa onko koordinaatit kunnollisia.
  % modified by Lu Cheng, 05.12.2012

  viallinen = 1;
  if ~isnumeric(coordinates)
    warning('Coordinates are not numerical!');
    return;
  end

  oikeanKokoinen = (size(coordinates,1) == ninds) & (size(coordinates,2) == 2);
  if ~oikeanKokoinen
    warning('Wrong coordinates dimension!');
    return;
  end

  posstr = cellfun(@(x) sprintf('%.10f',x),num2cell(coordinates),'UniformOutput',false);
  posstr = cellfun(@(x) regexprep(x,'0+$',''),posstr,'UniformOutput',false);

  uni1 = unique(posstr(:,1));
  uni2 = unique(posstr(:,2));
  posstr_new = posstr;

  if length(uni1)==ninds && length(uni2)==ninds
    viallinen = 0;
    return;
  else
    ans = questdlg('Input coordinates are not unique. Do you want to make them unique?','coordinates NOT unique', 'Yes','No','Yes');
    if strcmp(ans,'No')
      warning('Coordinates are not unique!');
      return;
    end
  end

  for i=1:length(uni1)
    tmpinds = find(ismember(posstr(:,1),uni1(i)));
    tmpNinds = length(tmpinds);

    if tmpNinds==1
      continue;
    end

    assert(tmpNinds<100);
    tmparr = round(linspace(0,99,tmpNinds));
    tmparr = tmparr(randperm(tmpNinds));

    for j=1:tmpNinds
      posstr_new{tmpinds(j),1}=sprintf('%s%02d',posstr{tmpinds(j),1},tmparr(j));
    end
  end

  for i=1:length(uni2)
    tmpinds = find(ismember(posstr(:,2),uni2(i)));
    tmpNinds = length(tmpinds);

    if tmpNinds==1
      continue;
    end

    assert(tmpNinds<100);
    tmparr = round(linspace(0,99,tmpNinds));
    tmparr = tmparr(randperm(tmpNinds));

    for j=1:tmpNinds
      posstr_new{tmpinds(j),2}=sprintf('%s%02d',posstr{tmpinds(j),2},tmparr(j));
    end
  end

  coordinates = cellfun(@str2double,posstr_new);
  uni1 = unique(coordinates(:,1));
  uni2 = unique(coordinates(:,2));
  if length(uni1)==ninds && length(uni2)==ninds
    viallinen = 0;
  else
    warning('Can not make coordinates unique!');
  end
end

function [cliques, separators, vorPoints, vorCells, pointers] ...
  = handleCoords(coordinates)
  %Laskee yksilצiden luonnolliset naapurit koordinaateista.
  %Naapurit lasketaan lisההmהll?koordinaatteihin pisteit?
  %jotta kutakin yksilצה vastaisi rajoitettu voronoi-solu
  %Puuttuvat koordinaatit (negatiiviset) tulevat erakkopisteiksi
  %
  %Mההrittהה lisהksi yksilצit?vastaavat voronoi tesselaation solut.
  %vorPoints:ssa on solujen kulmapisteet ja vorCells:ss?kunkin solun
  %kulmapisteiden indeksit. Pointers{i} sisהltהה solussa i olevien yksilצiden
  %indeksit.



  ninds = length(coordinates);
  [I,J] = find(coordinates>0 | coordinates <0);  %Kהsitellההn vain yksilצit? joilta koordinaatit
  I = unique(I);                %olemassa
  ncoords = length(I);
  puuttuvat = setdiff(1:ninds, I);
  new_coordinates = addPoints(coordinates(I,:)); %Ympהrצidההn yksilצt apupisteill?


  apuData = [new_coordinates(1:ncoords,:) (1:ncoords)'];
  apuData = sortrows(apuData,[1 2]);
  erot = [diff(apuData(:,1)) diff(apuData(:,2))];
  empties = find(erot(:,1)==0 & erot(:,2)==0);
  samat = cell(length(empties),1);
  pointer = 0;

  for i = 1:length(empties)
    if i == 1 | empties(i) - empties(i-1) > 1  %Tutkitaan onko eri pisteess?kuin edellinen
      pointer = pointer+1;
      samat{pointer} = [apuData(empties(i),3) apuData(empties(i)+1,3)];
    else
      samat{pointer} = [samat{pointer} apuData(empties(i)+1,3)];
    end
  end

  samat = samat(1:pointer);

  erot = []; apuData = []; empties = [];

  %tri = delaunay(new_coordinates(:,1), new_coordinates(:,2), {'Qt','Qbb','Qc','Qz'});    %Apupisteiden takia ok.
  tri = delaunay(new_coordinates(:,1), new_coordinates(:,2));
  %[rivi,sarake] = find(tri>ncoords);    %Jהtetההn huomiotta apupisteet
  %tri(rivi,:) = [];
  pituus = tri(:,1);
  pituus = length(pituus);
  parit = zeros(6*pituus,2);
  for i = 1:pituus                        %Muodostetaan kolmikoista parit
    j = 6*(i-1)+1;
    parit(j,:) = tri(i,1:2);
    parit(j+1,:) = tri(i,1:2:3);
    parit(j+2,:) = tri(i,2:3);
    parit(j+3:j+5,:) = [parit(j:j+2,2) parit(j:j+2,1)];
  end
  parit = unique(parit,'rows');
  [rivi,sarake] = find(parit>ncoords);     %Jהtetההn huomiotta apupisteet
  parit(rivi,:) = [];
  parit = I(parit);                         %Otetaan poistetut takaisin mukaan
  graph = sparse(parit(:,1),parit(:,2),1, ninds, ninds);


  %Kopioidaan samassa pisteess?olevien yksilצiden naapurustot
  %silt? jolle ne laitettu.

  for i = 1:length(samat);
    taulu = I(samat{i});
    [rivi,sarake] = find(graph(taulu,:)>0);
    if length(rivi) > 0
      kopioitava = graph(taulu(rivi(1)),:);
      for j = 1:length(taulu);
        graph(taulu(j),:) = kopioitava;
        graph(:,taulu(j)) = kopioitava';
      end
    end
  end

  %Asetetaan samassa pisteess?olevat yksilצt toistensa naapureiksi

  for i = 1:length(samat)
    for j = I(samat{i})
      for k = I(samat{i})
        if k ~= j
          graph(j,k) = 1;
        end
      end
    end
  end

  %Laskee maksimin klikkien ja separaattorien koolle
  %Mההritetההn myצs klikit ja separaattorit

  [ncliq, nsep, cliq, sep] = laskeKlikit(graph, ninds, ninds);

  sumcliq = sum(ncliq);
  sumsep = sum(nsep);
  maxCliqSize = max(find(sumcliq > 0));
  maxSepSize = max(find(sumsep > 0));

  cliques = zeros(length(cliq), maxCliqSize);
  separators = zeros(length(sep), maxSepSize);

  nollia = zeros(1, length(cliq));
  for i = 1:length(cliq);
    klikki = cliq{i};
    if length(klikki)>1
      cliques(i, 1:length(klikki)) = klikki;
    else
      nollia(i)=1;
    end
  end
  cliques(find(nollia==1), :) = [];

  for i = 1:length(sep);
    klikki = sep{i};
    separators(i, 1:length(klikki)) = klikki;
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %Mההritetההn yksilצit?vastaavat voronoi tesselaation solut

  [vorPoints, vorCells] = voronoin(new_coordinates, {'Qbb', 'Qz'});

  bounded = ones(length(vorCells),1);
  for i=1:length(vorCells)
    if isempty(vorCells{i}) || length(find(vorCells{i}==1))>0
      bounded(i)=0;
    end
  end



  vorCells = vorCells(find(bounded == 1));
  pointers = cell(length(vorCells),1);
  empties = zeros(1,length(vorCells));
  X = coordinates(:,1);
  Y = coordinates(:,2);

  for i=1:length(pointers)
    vx = vorPoints(vorCells{i},1);
    vy = vorPoints(vorCells{i},2);
    IN = inpolygon(X,Y,vx,vy);
    if any(IN)==0
      empties(i) = 1;
    else
      pointers{i} = find(IN ==1)';
    end
  end

  vorCells = vorCells(find(empties == 0));
  pointers = pointers(find(empties == 0));
end

function [ncliques, nseparators, cliques, separators] = laskeKlikit(M, maxCliqSize,maxSepSize)
  %Laskee samankokoisten klikkien mההrהn verkosta M
  %ncliques(i)=kokoa i olevien klikkien mההr?
  %nseparators vastaavasti

  ncliques=zeros(1,maxCliqSize);
  nseparators=zeros(1,maxSepSize);

  if isequal(M,[])
      return;
  end

  [cliques,separators]=findCliques(M);

  for i=1:length(cliques)
      ncliques(length(cliques{i}))=ncliques(length(cliques{i}))+1;
  end

  %cliqmax=max(find(ncliques~=0));
  %ncliques=ncliques(1:cliqmax);

  for i=1:length(separators)
      nseparators(length(separators{i}))=nseparators(length(separators{i}))+1;
  end
end
