
function kunnossa = testaaGenePopData(tiedostonNimi)
  % kunnossa == 0, jos data ei ole kelvollinen genePop data.
  % Muussa tapauksessa kunnossa == 1.

  kunnossa = 0;
  fid = fopen(tiedostonNimi);
  line1 = fgetl(fid);  %ensimmהinen rivi
  line2 = fgetl(fid);  %toinen rivi
  line3 = fgetl(fid);  %kolmas

  if (isequal(line1,-1) | isequal(line2,-1) | isequal(line3,-1))
    disp('Incorrect file format 1168'); fclose(fid);
    return
  end
  if (testaaPop(line1)==1 | testaaPop(line2)==1)
    disp('Incorrect file format 1172'); fclose(fid);
    return
  end
  if testaaPop(line3)==1
    %2 rivi tהllצin lokusrivi
    nloci = rivinSisaltamienMjonojenLkm(line2);
    line4 = fgetl(fid);
    if isequal(line4,-1)
      disp('Incorrect file format 1180'); fclose(fid);
      return
    end
    if ~any(line4==',')
      % Rivin nelj?tהytyy sisהltהה pilkku.
      disp('Incorrect file format 1185'); fclose(fid);
      return
    end
    pointer = 1;
    while ~isequal(line4(pointer),',')  %Tiedetההn, ett?pysהhtyy
      pointer = pointer+1;
    end
    line4 = line4(pointer+1:end);  %pilkun jהlkeinen osa
    nloci2 = rivinSisaltamienMjonojenLkm(line4);
    if (nloci2~=nloci)
      disp('Incorrect file format 1195'); fclose(fid);
      return
    end
  else
    line = fgetl(fid);
    lineNumb = 4;
    while (testaaPop(line)~=1 & ~isequal(line,-1))
      line = fgetl(fid);
      lineNumb = lineNumb+1;
    end
    if isequal(line,-1)
      disp('Incorrect file format 1206'); fclose(fid);
      return
    end
    nloci = lineNumb-2;
    line4 = fgetl(fid);  %Eka rivi pop sanan jהlkeen
    if isequal(line4,-1)
      disp('Incorrect file format 1212'); fclose(fid);
      return
    end
    if ~any(line4==',')
      % Rivin tהytyy sisהltהה pilkku.
      disp('Incorrect file format 1217'); fclose(fid);
      return
    end
    pointer = 1;
    while ~isequal(line4(pointer),',')  %Tiedetההn, ett?pysהhtyy.
      pointer = pointer+1;
    end

    line4 = line4(pointer+1:end);  %pilkun jהlkeinen osa
    nloci2 = rivinSisaltamienMjonojenLkm(line4);
    if (nloci2~=nloci)
      disp('Incorrect file format 1228'); fclose(fid);
      return
    end
  end
  kunnossa = 1;
  fclose(fid);
end
