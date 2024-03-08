function processed_data = process_data(filename, file_type, partitionCompare)
  switch file_type
    case 'BAPS'
      processed_data = process_BAPS_data(filename, partitionCompare);
    case 'GenePop'
      waitALittle;
      [filename, pathname] = uigetfile('*.txt', 'Load data in GenePop-format');
      if filename==0
        return;
      end
      if ~isempty(partitionCompare)
        fprintf(1,'Data: %s\n',[pathname filename]);
      end
      kunnossa = testaaGenePopData([pathname filename]);
      if kunnossa==0
        return
      end
      [data,popnames]=lueGenePopData([pathname filename]);

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
    case 'Preprocessed'
      waitALittle;
      [filename, pathname] = uigetfile('*.mat', 'Load pre-processed data');
      if filename==0
        return;
      end
      h0 = findobj('Tag','filename1_text');
      set(h0,'String',filename); clear h0;
      if ~isempty(partitionCompare)
        fprintf(1,'Data: %s\n',[pathname filename]);
      end

      struct_array = load([pathname filename]);
      if isfield(struct_array,'c')  %Matlab versio
        c = struct_array.c;
        if ~isfield(c,'dist')
          disp('Incorrect file format');
          return
        end
      elseif isfield(struct_array,'dist')  %Mideva versio
        c = struct_array;
      else
        disp('Incorrect file format');
        return;
      end
      processed_data = c;
    otherwise
      disp('Unknown file type');
  end
end
