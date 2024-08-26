
function [data, popnames] = lueGenePopData(tiedostonNimi)

  fid = fopen(tiedostonNimi);
  line = fgetl(fid);  %ensimmהinen rivi
  line = fgetl(fid);  %toinen rivi
  count = rivinSisaltamienMjonojenLkm(line);

  line = fgetl(fid);
  lokusRiveja = 1;
  while (testaaPop(line)==0)
    lokusRiveja = lokusRiveja+1;
    line = fgetl(fid);
  end

  if lokusRiveja>1
    nloci = lokusRiveja;
  else
    nloci = count;
  end

  popnames = cell(10,2);
  data = zeros(100, nloci+1);
  nimienLkm=0;
  ninds=0;
  poimiNimi=1;
  digitFormat = -1;
  while line ~= -1
    line = fgetl(fid);

    if poimiNimi==1
      %Edellinen rivi oli 'pop'
      nimienLkm = nimienLkm+1;
      ninds = ninds+1;
      if nimienLkm>size(popnames,1);
        popnames = [popnames; cell(10,2)];
      end
      nimi = lueNimi(line);
      if digitFormat == -1
        digitFormat = selvitaDigitFormat(line);
        divider = 10^digitFormat;
      end
      popnames{nimienLkm, 1} = {nimi};   %Nהin se on greedyMix:issהkin?!?
      popnames{nimienLkm, 2} = ninds;
      poimiNimi=0;

      data = addAlleles(data, ninds, line, divider);

    elseif testaaPop(line)
      poimiNimi = 1;

    elseif line ~= -1
      ninds = ninds+1;
      data = addAlleles(data, ninds, line, divider);
    end
  end

  data = data(1:ninds*2,:);
  popnames = popnames(1:nimienLkm,:);
  fclose(fid);
end
