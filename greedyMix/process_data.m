process_data(filename, file_type)
  switch file_type
    case 'BAPS-format'
    [filename, pathname] = uigetfile('*.txt', 'Load data in BAPS-format');
    if filename==0
        return;
    end

    if ~isempty(partitionCompare)
        fprintf(1,'Data: %s\n',[pathname filename]);
    end

    data = load([pathname filename]);
    ninds = testaaOnkoKunnollinenBapsData(data);  %TESTAUS
    if (ninds==0)
        disp('Incorrect Data-file.');
        return;
    end
    h0 = findobj('Tag','filename1_text');
    set(h0,'String',filename); clear h0;

    input_pops = questdlg(['When using data which are in BAPS-format, '...
        'you can specify the sampling populations of the individuals by '...
        'giving two additional files: one containing the names of the '...
        'populations, the other containing the indices of the first '...
        'individuals of the populations. Do you wish to specify the '...
        'sampling populations?'], ...
        'Specify sampling populations?',...
        'Yes', 'No', 'No');
    if isequal(input_pops,'Yes')
        waitALittle;
        [namefile, namepath] = uigetfile('*.txt', 'Load population names');
        if namefile==0
            kysyToinen = 0;
        else
            kysyToinen = 1;
        end
        if kysyToinen==1
            waitALittle;
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
%                 save(kokonimi,'c');
        save(kokonimi,'c','-v7.3'); % added by Lu Cheng, 08.06.2012
        clear c;
    end;

case 'GenePop-format'
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
%                 save(kokonimi,'c');
        save(kokonimi,'c','-v7.3'); % added by Lu Cheng, 08.06.2012
        clear c;
    end;

case 'Preprocessed data'
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
    data = double(c.data); rowsFromInd = c.rowsFromInd; alleleCodes = c.alleleCodes;
    noalle = c.noalle; adjprior = c.adjprior; priorTerm = c.priorTerm;
    dist = c.dist; popnames = c.popnames; Z = c.Z;
    clear c;
    otherwise
    return
  end
end
