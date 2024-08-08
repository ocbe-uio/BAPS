function processed_data = process_FASTA_data(file, partitionCompare)
  [filename1, pathname1] = uigetfile({'*.fasta';'*.*'}, 'Load data in FASTA-format');
  if filename1==0
    return;
  end

  if ~isempty(partitionCompare)
    fprintf(1,'Data: %s\n',[pathname filename]);
  end
  [heds, seqs] = fastaread([pathname1 filename1]);
  seqs = seqs(:);
  alnMat = cell2mat(seqs);
  nSeq = length(seqs);
  clear seqs;

  cc = preprocAln(alnMat);
  cc.heds = heds;
  cc.nSeq = nSeq;
  dist = seqpdist(alnMat,'method','p-distance');
  Z = linkage(dist,'complete');

  ninds = nSeq;

  clear alnMat heds nSeq

  [filename2,pathname2]=uigetfile('*.txt', 'Load individual coordinates');
  if filename2==0
    return
  end

  coordinates = load([pathname2 filename2]);
  [viallinen coordinates] = testaaKoordinaatit(ninds, coordinates); % added by Lu Cheng, 05.12.2012
  if viallinen
    disp('Incorrect coordinates');
    return
  end

  inp = [filename1 ' & ' filename2];
  h0 = findobj('Tag','filename1_text');
  set(h0,'String',inp);
  clear h0; clear inp;
  clear filename1; clear filename2; clear pathname1; clear pathname2;

  input_pops = questdlg(['When using data which are in FASTA-format, '...
  'you can specify the sampling populations of the individuals by '...
  'giving two additional files: one containing the names of the '...
  'populations, the other containing the indices of the first '...
  'individuals of the populations. Do you wish to specify the '...
  'sampling populations?'], ...
  'Specify sampling populations?',...
  'Yes', 'No', 'No');
  if isequal(input_pops,'Yes')
    [namefile, namepath] = uigetfile('*.txt', 'Load population names');
    if namefile==0
      kysyToinen = 0;
    else
      kysyToinen = 1;
    end
    if kysyToinen==1
      [indicesfile, indicespath] = uigetfile('*.txt', 'Load population indices');
      if indicesfile==0
        popnames = [];
      else
        popnames = initPopNames([namepath namefile],[indicespath indicesfile]);
      end
    else
      popnames = [];
    end
  else
    popnames = [];
  end

  disp('Pre-processing the data. This may take several minutes.');

  [cliques, separators, vorPoints, vorCells, pointers] = ...
  handleCoords(coordinates);

  cc.locCliques = cliques;
  cc.locSeparators = separators;
  cc.popnames = popnames;
  cc.vorPoints = vorPoints;
  cc.vorCells = vorCells;
  cc.pointers = pointers;
  cc.coordinates = coordinates;
  format_type = 'FASTA';

  save_preproc = questdlg('Do you wish to save pre-processed data?',...
  'Save pre-processed data?',...
  'Yes','No','Yes');
  if isequal(save_preproc,'Yes')
    [filename, pathname] = uiputfile('*.mat','Save pre-processed data as');
    kokonimi = [pathname filename];
    save(kokonimi,'cc','dist','Z','format_type','-v7.3'); % added by Lu Cheng, 08.06.2012
  end

  handleIndiFastaCase(cc,dist,Z);

end
