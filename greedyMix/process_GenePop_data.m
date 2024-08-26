
function processed_data = process_GenePop_data(filename)
  kunnossa = testaaGenePopData(filename);
  if kunnossa==0
    return
  end
  [data,popnames]=lueGenePopData(filename);

  h0 = findobj('Tag','filename1_text');
  set(h0,'String',filename); clear h0;

  [data, rowsFromInd, alleleCodes, noalle, adjprior, priorTerm] = handleData(data);
  [Z,dist] = newGetDistances(data,rowsFromInd);

  % Forming and returning pre-processed data
  c.data = data; c.rowsFromInd = rowsFromInd; c.alleleCodes = alleleCodes;
  c.noalle = noalle; c.adjprior = adjprior; c.priorTerm = priorTerm;
  c.dist = dist; c.popnames = popnames; c.Z = Z;
  processed_data = c;
end
