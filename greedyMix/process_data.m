function processed_data = process_data(filename, file_type, partitionCompare, coordinates)
  switch file_type
    case 'BAPS'
      if filename == ""
        [filename, pathname] = uigetfile('*.txt', 'Load data in BAPS-format');
        filename = fullfile(pathname, filename);
      end
      processed_data = process_BAPS_data(filename, partitionCompare);
    case 'FASTA'
      processed_data = process_FASTA_data(filename, partitionCompare, coordinates);
    case 'GenePop'
      if filename == ""
        [filename, pathname] = uigetfile('*.txt', 'Load data in GenePop-format');
        filename = fullfile(pathname, filename);
      end
      if ~isempty(partitionCompare)
        fprintf(1,'Data: %s\n', filename);
      end
      processed_data = process_GenePop_data(filename);
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
      fprintf('Unknown file type: %s\n', file_type);
  end
end
