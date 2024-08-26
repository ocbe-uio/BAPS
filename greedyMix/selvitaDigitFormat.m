function df = selvitaDigitFormat(line)
  % line on ensimmהinen pop-sanan jהlkeinen rivi
  % Genepop-formaatissa olevasta datasta. funktio selvittהה
  % rivin muodon perusteella, ovatko datan alleelit annettu
  % 2 vai 3 numeron avulla.

  n = 1;
  merkki = line(n);
  while ~isequal(merkki,',')
    n = n+1;
    merkki = line(n);
  end

  while ~any(merkki == '0123456789');
    n = n+1;
    merkki = line(n);
  end
  numeroja = 0;
  while any(merkki == '0123456789');
    numeroja = numeroja+1;
    n = n+1;
    merkki = line(n);
  end

  df = numeroja/2;
end
