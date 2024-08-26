function data = addAlleles(data, ind, line, divider)
  % Lisaa BAPS-formaatissa olevaan datataulukkoon
  % yksilצה ind vastaavat rivit. Yksilצn alleelit
  % luetaan genepop-formaatissa olevasta rivist?
  % line. Jos data on 3 digit formaatissa on divider=1000.
  % Jos data on 2 digit formaatissa on divider=100.

  nloci = size(data,2)-1;
  if size(data,1) < 2*ind
    data = [data; zeros(100,nloci+1)];
  end

  k=1;
  merkki=line(k);
  while ~isequal(merkki,',')
    k=k+1;
    merkki=line(k);
  end
  line = line(k+1:end);
  clear k; clear merkki;

  alleeliTaulu = sscanf(line,'%d');

  if length(alleeliTaulu)~=nloci
    disp('Incorrect data format.');
  end

  for j=1:nloci
    ekaAlleeli = floor(alleeliTaulu(j)/divider);
    if ekaAlleeli==0;
      ekaAlleeli=-999;
    end
    tokaAlleeli = rem(alleeliTaulu(j),divider);
    if tokaAlleeli==0;
      tokaAlleeli=-999;
    end

    data(2*ind-1,j) = ekaAlleeli;
    data(2*ind,j) = tokaAlleeli;
  end

  data(2*ind-1,end) = ind;
  data(2*ind,end) = ind;
end
