
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
  save_preproc = questdlg('Do you wish to save pre-processed data?',...
  'Save pre-processed data?',...
  'Yes','No','Yes');
  if isequal(save_preproc,'Yes');
    waitALittle;
    [filename, pathname] = uiputfile('*.mat','Save pre-processed data as');
    kokonimi = [pathname filename];
    c.data = data; c.rowsFromInd = rowsFromInd; c.alleleCodes = alleleCodes;
    c.noalle = noalle; c.adjprior = adjprior; c.priorTerm = priorTerm;
    c.dist = dist; c.popnames = popnames; c.Z = Z;
    save(kokonimi,'c','-v7.3'); % added by Lu Cheng, 08.06.2012
    processed_data = c;
  end
