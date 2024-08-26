function nimi = lueNimi(line)
  %Palauttaa line:n alusta sen osan, joka on ennen pilkkua.
  n = 1;
  merkki = line(n);
  nimi = '';
  while ~isequal(merkki,',')
    nimi = [nimi merkki];
    n = n+1;
    merkki = line(n);
  end
end
