function diffInCounts = computeDiffInCounts(rows, max_noalle, nloci, data)
  % Muodostaa max_noalle*nloci taulukon, jossa on niiden alleelien
  % lukumההrהt (vastaavasti kuin COUNTS:issa), jotka ovat data:n
  % riveillה rows. rows pitהה olla vaakavektori.

  diffInCounts = zeros(max_noalle, nloci);
  for i=rows
    row = data(i, :);
    notEmpty = find(row >= 0);

    if ~isempty(notEmpty)
      element = row(notEmpty) + (notEmpty-1) * max_noalle;
      diffInCounts(element) = diffInCounts(element) + 1;
    end
  end
end
