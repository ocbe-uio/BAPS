function count = rivinSisaltamienMjonojenLkm(line)
  % Palauttaa line:n sisהltהmien mjonojen lukumההrהn.
  % Mjonojen vהliss?tהytyy olla vהlilyצnti.
  count = 0;
  pit = length(line);
  tila = 0;    %0, jos odotetaan vהlilyצntej? 1 jos odotetaan muita merkkej?
  for i=1:pit
    merkki = line(i);
    if (isspace(merkki) & tila==0)
      %Ei tehd?mitההn.
    elseif (isspace(merkki) & tila==1)
      tila = 0;
    elseif (~isspace(merkki) & tila==0)
      tila = 1;
      count = count+1;
    elseif (~isspace(merkki) & tila==1)
      %Ei tehd?mitההn
    end
  end
end
