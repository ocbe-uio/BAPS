function ninds = testaaOnkoKunnollinenBapsData(data)
  %Tarkastaa onko viimeisess?sarakkeessa kaikki
  %luvut 1,2,...,n johonkin n:ההn asti.
  %Tarkastaa lisהksi, ett?on vהhintההn 2 saraketta.
  if size(data, 1) < 2
    ninds = 0;
    return;
  end
  lastCol = data(:, end);
  ninds = max(lastCol);
  if ~isequal((1:ninds)', unique(lastCol))
    ninds = 0;
    return;
  end
end
