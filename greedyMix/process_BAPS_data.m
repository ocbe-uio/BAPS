function processed_data = process_BAPS_data(file, partitionCompare)
  if ~isempty(partitionCompare)
    fprintf(1, 'Data: %s\n', file);
  end
  data = importdata(file);
  ninds = testaaOnkoKunnollinenBapsData(data);  % for testing purposes?
  if (ninds == 0)
    warning('Incorrect Data-file.');
    return;
  end
  popnames = []; % Dropped specification of population names (from BAPS 6)

  [data, rowsFromInd, alleleCodes, noalle, adjprior, priorTerm] = handleData(data);
  [Z, dist] = newGetDistances(data, rowsFromInd);

  % Forming and saving pre-processed data
  processed_data.data = data;
  processed_data.rowsFromInd = rowsFromInd;
  processed_data.alleleCodes = alleleCodes;
  processed_data.noalle = noalle;
  processed_data.adjprior = adjprior;
  processed_data.priorTerm = priorTerm;
  processed_data.dist = dist;
  processed_data.popnames = popnames;
  processed_data.Z = Z;
end
