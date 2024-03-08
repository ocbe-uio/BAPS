function c_result = semi_linkageMix(c_train, c_test, npopsTable)
% Greedy search algorithm with unknown number of classes for linkage
% clustering.

% Modified from linkageMix.m by Lu Cheng, 16.02.2010
% Input: preprocessed training data and test data, as well as the
% population number, npopsTable is a vector
% NOTE, here we assume that c_train contains all the correct priors

% IS_SINGLE_SAMPLE is introduced to handle the case of only one test sample
% Lu Cheng, 07.03.2011

global POP_LOGML;       global PARTITION;
global CQ_COUNTS;       global SP_COUNTS;       %These counts are for populations
global CQ_SUMCOUNTS;    global SP_SUMCOUNTS;    %not for individuals
global LOGDIFF;

global TRAIN_CQ_COUNTS;
global TRAIN_SP_COUNTS;
global TRAIN_CQ_SUMCOUNTS;
global TRAIN_SP_SUMCOUNTS;

clearGlobalVars;

% noalle = c_train.noalle;
% adjprior = c_train.adjprior; %priorTerm = c.priorTerm;
% rowsFromInd = c_train.rowsFromInd;
adjprior_cq = c_train.adjprior_cq;
adjprior_sp = c_train.adjprior_sp;

counts_cq = c_test.counts_cq;
counts_sp = c_test.counts_sp;

% UPDATE the priors by the following constants from the training data
% Lu Cheng, 17.02.2010 ------------------------ %
clusters = unique(c_train.cluster_labels);
n_train_clusters = length(clusters);
if ~all(clusters'==1:n_train_clusters) %labels should cover 1 to n_train_clusters
    error(cat(2,'Error! Cluster labels should range from 1 to ',n_train_clusters,'!'));
end
[d1_cq d2_cq d3_cq] = size(c_train.counts_cq);
[d1_sp d2_sp d3_sp] = size(c_train.counts_sp);
TRAIN_CQ_COUNTS = zeros(d1_cq,d2_cq,n_train_clusters);
TRAIN_SP_COUNTS = zeros(d1_sp,d2_sp,n_train_clusters);
for i = 1:n_train_clusters
    cluster_inds = (c_train.cluster_labels == clusters(i));
    TRAIN_CQ_COUNTS(:,:,i) = sum(c_train.counts_cq(:,:,cluster_inds),3);
    TRAIN_SP_COUNTS(:,:,i) = sum(c_train.counts_sp(:,:,cluster_inds),3);
end

train_counts_cq = TRAIN_CQ_COUNTS;
train_counts_sp = TRAIN_SP_COUNTS;

clear clusters cluster_inds i
clear d1_cq d2_cq d3_cq d1_sp d2_sp d3_sp cluster_inds
% --------------------------------------------- %

IS_SINGLE_SAMPLE = false;
if isfield(c_test,'dist')
    dist = c_test.dist; Z = c_test.Z;
else
    IS_SINGLE_SAMPLE = true;
end
clear c_test c_train;

nruns = length(npopsTable);
logmlBest = -1e50;
partitionSummary = -1e50*ones(100,2);  % 100 best partitions (npops and logml)
partitionSummary(:,1) = zeros(100,1);
worstLogml = -1e50; worstIndex = 1;

for run = 1:nruns
    if IS_SINGLE_SAMPLE && (npopsTable(run)>n_train_clusters+1)
        npops = n_train_clusters + 1;
    else
        npops = npopsTable(run);
    end

    dispLine;
    disp(['Run ' num2str(run) '/' num2str(nruns) ...
        ', maximum number of populations ' num2str(npops) '.']);

    %----- added by Lu Cheng, 23.02.2010------%
    % modify the training data so it fits the dimension of given population
    % number
    %n_train_clusters = size(train_counts_cq,3);
    if npops<n_train_clusters
        TRAIN_CQ_COUNTS = train_counts_cq(:,:,1:npops);
        TRAIN_SP_COUNTS = train_counts_sp(:,:,1:npops);
    elseif npops==n_train_clusters
        TRAIN_CQ_COUNTS = train_counts_cq;
        TRAIN_SP_COUNTS = train_counts_sp;
    else
        [r1 c1 p1] = size(train_counts_cq);
        [r2 c2 p2] = size(train_counts_sp);
        TRAIN_CQ_COUNTS = zeros(r1, c1, npops);
        TRAIN_SP_COUNTS = zeros(r2, c2, npops);
        TRAIN_CQ_COUNTS(:,:,1:n_train_clusters) = train_counts_cq;
        TRAIN_SP_COUNTS(:,:,1:n_train_clusters) = train_counts_sp;
    end

    TRAIN_CQ_SUMCOUNTS = squeeze(sum(TRAIN_CQ_COUNTS(:,:,:),1))';
    TRAIN_SP_SUMCOUNTS = squeeze(sum(TRAIN_SP_COUNTS(:,:,:),1))';

    clear r1 c1 p1 r2 c2 p2
    %---------------------------------------%

    if IS_SINGLE_SAMPLE
        PARTITION = 1;
    else
        PARTITION = admixture_initialization_2(npops, Z);
    end

    [cq_counts, cq_sumcounts] = calCounts(counts_cq, npops);
    CQ_COUNTS = cq_counts;  clear cq_counts;
    CQ_SUMCOUNTS = cq_sumcounts; clear cq_sumcounts;

    [sp_counts, sp_sumcounts] = calCounts(counts_sp,npops);
    SP_COUNTS = sp_counts;  clear sp_counts;
    SP_SUMCOUNTS = sp_sumcounts; clear sp_sumcounts;

    logml = computeLogml(adjprior_cq, adjprior_sp);
    POP_LOGML = computePopulationLogml(1:npops,adjprior_cq, adjprior_sp);

    %disp(POP_LOGML');

    if logml>worstLogml
        [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex, n_train_clusters);
        if (added==1)
            [worstLogml, worstIndex] = min(partitionSummary(:,2));
        end
    end

    %%%%
    % Finding the best partition with the greedy search algorithm
    %%%%
    nRoundTypes = 7;
    tested = zeros(nRoundTypes,1);
    roundTypes = [1 1];
    ready = 0;  phase = 1;
    ninds = length(PARTITION);  % number of individuals
    LOGDIFF = repmat(-Inf,ninds,npops);

    disp(' ');
    disp(['Mixture analysis started with initial ' num2str(npops) ' populations.']);

    while ready ~= 1
        changesMade = 0;

        disp(['Performing steps: ' num2str(roundTypes)]);

        for n = 1:length(roundTypes)
            round = roundTypes(n);
            % pack;
            if tested(round) == 1

            elseif round==1   % Moving one individual to another population
                inds = randperm(ninds);     % random order
                changesMadeNow = 0;
                for ind = inds
                    i1 = PARTITION(ind);
                    indCqCounts = uint16(counts_cq(:,:,ind));
                    indSpCounts = uint16(counts_sp(:,:,ind));
                    changesInLogml = computeChanges(ind, adjprior_cq, ...
                        adjprior_sp, indCqCounts, indSpCounts);

                    [maxChange, i2] = max(changesInLogml);

                    if (i1~=i2 && maxChange>1e-5)

                        %disp(changesInLogml);

                        % Individual is moved
                        changesMade = 1;
                        if changesMadeNow == 0
                            disp('action 1');
                            changesMadeNow = 1;
                            tested = zeros(nRoundTypes,1);
                        end
                        updateGlobalVariables(ind, i2, indCqCounts, ...
                            indSpCounts, adjprior_cq, adjprior_sp);
                        logml = logml+maxChange;

                        if logml>worstLogml
                            [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex, n_train_clusters);
                            if (added==1)
                                [worstLogml, worstIndex] = min(partitionSummary(:,2));
                            end
                        end
                    end
                end

                if changesMadeNow == 0
                    tested(round) = 1;
                end


                %disp(PARTITION');  % for test

            elseif round==2  % Combining two populations
                maxChange = 0;

                %n_train_clusters % see line 32, Lu Cheng, 04.03.2010

%                 for pop = 1:npops
%                     changesInLogml = computeChanges2(pop, adjprior_cq, adjprior_sp); % all inds in 'pop' is moved to other clusters
%                     [biggest, index] = max(changesInLogml);
%                     if biggest>maxChange
%                         maxChange = biggest;
%                         i1 = pop;
%                         i2 = index;
%                         fprintf('moving population %d to population %d.\n',i1,i2);
%                     end
%                     disp(changesInLogml');
%                 end

                % modified by Lu Cheng, 04.03.2010
                % here we only combine the outer clusters with other
                % outerclusters or training clusters. In case there is only
                % one outer cluster in the PARTITION, it will be combined with
                % some training cluster, depending on the logml
                if npops>n_train_clusters
                    for pop = n_train_clusters+1:npops
                        if any(PARTITION==pop)
                            changesInLogml = computeChanges2(pop, adjprior_cq, adjprior_sp); % all inds in 'pop' is moved to other clusters
                            [biggest, index] = max(changesInLogml);
                            if biggest>maxChange
                                maxChange = biggest;
                                i1 = pop;
                                i2 = index;
                                %fprintf('moving population %d to population %d.\n',i1,i2);
                            end
                        end
                    end

                    %fprintf('maxChange: %d\n',maxChange);
                end


                if maxChange>1e-5
                    disp('action 2');
                    changesMade = 1;
                    tested = zeros(nRoundTypes,1);

                    updateGlobalVariables2(i1, i2, adjprior_cq, adjprior_sp); %all inds in i1 are moved to i2
                    logml = logml + maxChange;
                    if logml>worstLogml
                        [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex, n_train_clusters);
                        if (added==1)
                            [worstLogml, worstIndex] = min(partitionSummary(:,2));
                        end
                    end

                    %disp(PARTITION');  % for test

                else
                    tested(round) = 1;
                end


            elseif (round==3 || round==4) && ~IS_SINGLE_SAMPLE  % Splitting population into smaller groups
                maxChange = 0;

                for pop = 1:npops
                    inds2 = find(PARTITION==pop);
                    ninds2 = length(inds2);
                    if ninds2>5
                        % Computing the distance between individuals inds2
                        dist2 = laskeOsaDist(inds2, dist, ninds);
                        Z2 = computeLinkage(dist2');

                        % Number of groups:
                        if round==3
                            npops2 = min(20, floor(ninds2 / 5));
                        elseif round==4
                            npops2 = 2;
                        end
                        T2 = cluster_own(Z2, npops2);

                        changesInLogml = computeChanges3(T2, inds2, pop, ...
                            counts_cq, counts_sp, adjprior_cq, adjprior_sp); % npops2 * npops matrix
                        [biggest, index] = max(changesInLogml(1:end));
                        if biggest > maxChange
                            maxChange = biggest;
                            movingGroup = rem(index,npops2); % The group, which is moved
                            if movingGroup==0, movingGroup = npops2; end
                            movingInds = inds2(logical(T2==movingGroup));
                            i2 = ceil(index/npops2);    % pop where movingGroup would be moved
                        end
                    end
                end
                if maxChange>1e-5
                    changesMade = 1;
                    tested = zeros(nRoundTypes,1);
                    if round==3
                        disp('action 3');
                    else
                        disp('action 4');
                    end
                    indCqCounts = uint16(sum(counts_cq(:,:,movingInds),3));
                    indSpCounts = uint16(sum(counts_sp(:,:,movingInds),3));
                    updateGlobalVariables3(movingInds, i2,indCqCounts, ...
                        indSpCounts, adjprior_cq, adjprior_sp);
                    logml = logml + maxChange;

                    if logml>worstLogml
                        [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex, n_train_clusters);
                        if (added==1)
                            [worstLogml, worstIndex] = min(partitionSummary(:,2));
                        end
                    end
                else
                    tested(round) = 1;
                end

            elseif (round == 5 || round == 6) && ~IS_SINGLE_SAMPLE
                %Moving individuals out of population until positive change
                %in logml has occured
                pop=0;
                changesMadeNow = 0;
                %Saving old values
                poplogml = POP_LOGML;
                partition = PARTITION;
                cq_counts = CQ_COUNTS;
                sp_counts = SP_COUNTS;
                cq_sumcounts = CQ_SUMCOUNTS;
                sp_sumcounts = SP_SUMCOUNTS;
                logdiff = LOGDIFF;

                while (pop < npops && changesMadeNow == 0)
                    pop = pop+1;
                    totalChangeInLogml = 0;
                    inds = find(PARTITION==pop);
                    if round == 5
                        %Random order
                        aputaulu = [inds rand(length(inds),1)];
                        aputaulu = sortrows(aputaulu,2);
                        inds = aputaulu(:,1)';
                    elseif round == 6
                        inds = returnInOrder(inds, pop, counts_cq, counts_sp, ...
                            adjprior_cq, adjprior_sp);
                    end

                    i=0;

                    while (length(inds)>0 && i<length(inds))
                        i = i+1;
                        ind = inds(i);
                        indCqCounts = uint16(counts_cq(:,:,ind));
                        indSpCounts = uint16(counts_sp(:,:,ind));
                        changesInLogml = computeChanges(ind, adjprior_cq, ...
                            adjprior_sp, indCqCounts, indSpCounts);
                        changesInLogml(pop) = -1e50;   % Varmasti ei suurin!!!
                        [maxChange, i2] = max(changesInLogml);
                        updateGlobalVariables(ind, i2, indCqCounts, ...
                            indSpCounts, adjprior_cq, adjprior_sp);
                        totalChangeInLogml = totalChangeInLogml+maxChange;
                        logml = logml+maxChange;
                        if round == 6
                            % Stop immediatly when change in logml is
                            % positive
                            if totalChangeInLogml > 1e-5
                                i=length(inds);
                            end
                        end
                    end

                    if totalChangeInLogml>1e-5
                        if round == 5
                            disp('action 5');
                        elseif round == 6
                            disp('action 6');
                        end
                        tested = zeros(nRoundTypes,1);
                        changesMadeNow=1;
                        changesMade = 1;
                        if logml>worstLogml
                            [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex, n_train_clusters);
                            if (added==1)
                                [worstLogml, worstIndex] = min(partitionSummary(:,2));
                            end
                        end
                    else
                        % No better partition was found, restoring the old
                        % values
                        PARTITION = partition;
                        POP_LOGML = poplogml;
                        CQ_COUNTS = cq_counts;
                        SP_COUNTS = sp_counts;
                        CQ_SUMCOUNTS = cq_sumcounts;
                        SP_SUMCOUNTS = sp_sumcounts;
                        LOGDIFF = logdiff;
                        logml = logml - totalChangeInLogml;
                    end
                end
                clear partition;  clear poplogml;
                if changesMadeNow == 0
                    tested(round) = 1;
                end

            elseif round == 7 && ~IS_SINGLE_SAMPLE
                emptyPop = findEmptyPop(npops);
                j = 0;
                pops = randperm(npops);
                % totalChangeInLogml = 0;
                if emptyPop == -1
                    j = npops;
                end
                changesMadeNow = 0;
                while (j < npops)
                    j = j +1;
                    pop = pops(j);
                    inds2 = find(PARTITION == pop);
                    ninds2 = length(inds2);
                    if ninds2 > 5
                        partition = PARTITION;
                        cq_sumcounts = CQ_SUMCOUNTS;
                        cq_counts = CQ_COUNTS;
                        sp_sumcounts = SP_SUMCOUNTS;
                        sp_counts = SP_COUNTS;
                        poplogml = POP_LOGML;
                        logdiff = LOGDIFF;
                        % pack;

                        dist2 = laskeOsaDist(inds2, dist, ninds);
                        Z2 = computeLinkage(dist2');
                        T2 = cluster_own(Z2, 2);
                        % movingInds = inds2(find(T2 == 1));
                        movingInds = inds2(logical(T2 == 1));
                        changesInLogml = computeChanges3(T2, inds2, pop, ...
                            counts_cq, counts_sp, adjprior_cq, adjprior_sp);
                        totalChangeInLogml = changesInLogml(1, emptyPop);

                        indCqCounts = uint16(sum(counts_cq(:,:,movingInds),3));
                        indSpCounts = uint16(sum(counts_sp(:,:,movingInds),3));
                        updateGlobalVariables3(movingInds, emptyPop,indCqCounts, ...
                            indSpCounts, adjprior_cq, adjprior_sp);

                        changed = 1;

                        while (changed == 1)
                            changed = 0;

                            changesInLogml = computeChanges5(inds2, pop, emptyPop, ...
                                counts_cq, counts_sp, adjprior_cq, adjprior_sp);
                            [maxChange, index] = max(changesInLogml);
                            moving = inds2(index);
                            if (PARTITION(moving) == pop)
                                i2 = emptyPop;
                            else
                                i2 = pop;
                            end

                            if maxChange > 1e-5
                                indCqCounts = uint16(counts_cq(:,:,moving));
                                indSpCounts = uint16(counts_sp(:,:,moving));
                                updateGlobalVariables3(moving, i2,indCqCounts, ...
                                    indSpCounts, adjprior_cq, adjprior_sp);
                                changed = 1;
                                totalChangeInLogml = totalChangeInLogml + maxChange;
                            end
                        end

                        if totalChangeInLogml > 1e-5
                            changesMade = 1;
                            changesMadeNow = 1;
                            logml = logml + totalChangeInLogml;
                            if logml>worstLogml
                                [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex, n_train_clusters);
                                if (added==1)
                                    [worstLogml, worstIndex] = min(partitionSummary(:,2));
                                end
                            end
                            disp('action 7');
                            tested = zeros(nRoundTypes, 1);
                            j = npops;
                        else
                            % No better partition was found, restoring the old
                            % values
                            PARTITION = partition;
                            POP_LOGML = poplogml;
                            CQ_COUNTS = cq_counts;
                            SP_COUNTS = sp_counts;
                            CQ_SUMCOUNTS = cq_sumcounts;
                            SP_SUMCOUNTS = sp_sumcounts;
                            LOGDIFF = logdiff;
                            %logml = logml - totalChangeInLogml;
                        end
                    end
                end
                if changesMadeNow == 0
                    tested(round) = 1;
                end
            end

        end


        if changesMade == 0
            if phase==1
                phase = 2;
            elseif phase==2
                phase = 3;
            elseif phase==3
                phase = 4;
            elseif phase==4;
                phase = 5;
            elseif phase==5
                ready = 1;
            end
        else
            changesMade = 0;
        end

        if ready==0
            if phase==1
                roundTypes=[1];
            elseif phase==2
                roundTypes=[2];
            elseif phase==3
                roundTypes=[5 5 7];
            elseif phase==4
                roundTypes=[4 3 1 1];
            elseif phase==5
                roundTypes=[6 2 7 3 4 1];
            end
        end

    end
    % Saving results

    npops = removeEmptyPops(n_train_clusters,npops);
    POP_LOGML = computePopulationLogml(1:npops, adjprior_cq, adjprior_sp);

    disp(['Found partition with ' num2str(npops) ' populations.']);
    disp(['Log(ml) = ' num2str(logml)]);
    disp(' ');

    if logml>logmlBest
        % Updating the best found partition
        logmlBest = logml;
        npopsBest = npops;
        partitionBest = PARTITION;
        cq_countsBest = CQ_COUNTS;
        sp_countsBest = SP_COUNTS;
        cq_sumcountsBest = CQ_SUMCOUNTS;
        sp_sumcountsBest = SP_SUMCOUNTS;
        pop_logmlBest = POP_LOGML;
        logdiffbest = LOGDIFF;
    end
end

c_result.logml = logmlBest;
c_result.npops = npopsBest;
c_result.PARTITION = partitionBest;
c_result.CQ_COUNTS = cq_countsBest;
c_result.SP_COUNTS = sp_countsBest;
c_result.CQ_SUMCOUNTS = cq_sumcountsBest;
c_result.SP_SUMCOUNTS = sp_sumcountsBest;
c_result.POP_LOGML = pop_logmlBest;
c_result.LOGDIFF = logdiffbest;
c_result.partitionSummary = partitionSummary;

% Calculate the posterior probabilities if a sample i is moved to cluster
% j, from LOGDIFF (nsample*npops matrix)
% Each row of LOGDIFF represents logml changes if the sample is moved to
% the corresponding cluster
% See page 24 in BAPS manual
% added by Lu Cheng, 29.03.2010
[ninds npops] = size(c_result.LOGDIFF);
clusterProbTable = zeros(ninds, npops);
for i=1:ninds
    tmpRow = exp(c_result.LOGDIFF(i,:));
    clusterProbTable(i,:) = tmpRow/sum(tmpRow);
end
c_result.clusterProbTable = clusterProbTable;
% ------------

%--------------------------------------------------------------------------
% The next three functions are for computing the initial partition
% according to the distance between the individuals

function initial_partition = admixture_initialization_2(nclusters,Z)
initial_partition = cluster_own(Z,nclusters);

%--------------------------------------------------------------------------
function T = cluster_own(Z,nclust)
% true=logical(1);
% false=logical(0);

maxclust = nclust;
% Start of algorithm
m = size(Z,1)+1;
T = zeros(m,1);
% maximum number of clusters based on inconsistency
if m <= maxclust
    T = (1:m)';
elseif maxclust==1
    T = ones(m,1);
else
    clsnum = 1;
    for k = (m-maxclust+1):(m-1)
        i = Z(k,1); % left tree
        if i <= m % original node, no leafs
            T(i) = clsnum;
            clsnum = clsnum + 1;
        elseif i < (2*m-maxclust+1) % created before cutoff, search down the tree
            T = clusternum(Z, T, i-m, clsnum);
            clsnum = clsnum + 1;
        end
        i = Z(k,2); % right tree
        if i <= m  % original node, no leafs
            T(i) = clsnum;
            clsnum = clsnum + 1;
        elseif i < (2*m-maxclust+1) % created before cutoff, search down the tree
            T = clusternum(Z, T, i-m, clsnum);
            clsnum = clsnum + 1;
        end
    end
end

function T = clusternum(X, T, k, c)
m = size(X,1)+1;
while(~isempty(k))
    % Get the children of nodes at this level
    children = X(k,1:2);
    children = children(:);

    % Assign this node number to leaf children
    t = (children<=m);
    T(children(t)) = c;

    % Move to next level
    k = children(~t) - m;
end

%--------------------------------------------------------------------------

function dist2 = laskeOsaDist(inds2, dist, ninds)
% Muodostaa dist vektorista osavektorin, joka sis�lt�� yksil�iden inds2
% v�liset et�isyydet. ninds=kaikkien yksil�iden lukum��r�.

ninds2 = length(inds2);
apu = zeros(nchoosek(ninds2,2),2);
rivi = 1;
for i=1:ninds2-1
    for j=i+1:ninds2
        apu(rivi, 1) = inds2(i);
        apu(rivi, 2) = inds2(j);
        rivi = rivi+1;
    end
end
apu = (apu(:,1)-1).*ninds - apu(:,1) ./ 2 .* (apu(:,1)-1) + (apu(:,2)-apu(:,1));
dist2 = dist(apu);

%--------------------------------------------------------------------------

function Z = computeLinkage(Y, method)
[k, n] = size(Y);
m = (1+sqrt(1+8*n))/2;
if k ~= 1 || m ~= fix(m)
    error('The first input has to match the output of the PDIST function in size.');
end
if nargin == 1 % set default switch to be 'co'
    method = 'co';
end
method = lower(method(1:2)); % simplify the switch string.
% monotonic = 1;
Z = zeros(m-1,3); % allocate the output matrix.
N = zeros(1,2*m-1);
N(1:m) = 1;
n = m; % since m is changing, we need to save m in n.
R = 1:n;
for s = 1:(n-1)
    X = Y;
    [v, k] = min(X);
    i = floor(m+1/2-sqrt(m^2-m+1/4-2*(k-1)));
    j = k - (i-1)*(m-i/2)+i;
    Z(s,:) = [R(i) R(j) v]; % update one more row to the output matrix A
    I1 = 1:(i-1); I2 = (i+1):(j-1); I3 = (j+1):m; % these are temp variables.
    U = [I1 I2 I3];
    I = [I1.*(m-(I1+1)/2)-m+i i*(m-(i+1)/2)-m+I2 i*(m-(i+1)/2)-m+I3];
    J = [I1.*(m-(I1+1)/2)-m+j I2.*(m-(I2+1)/2)-m+j j*(m-(j+1)/2)-m+I3];

    switch method
        case 'si' %single linkage
            Y(I) = min(Y(I),Y(J));
        case 'av' % average linkage
            Y(I) = Y(I) + Y(J);
        case 'co' %complete linkage
            Y(I) = max(Y(I),Y(J));
        case 'ce' % centroid linkage
            K = N(R(i))+N(R(j));
            Y(I) = (N(R(i)).*Y(I)+N(R(j)).*Y(J)-(N(R(i)).*N(R(j))*v^2)./K)./K;
        case 'wa'
            Y(I) = ((N(R(U))+N(R(i))).*Y(I) + (N(R(U))+N(R(j))).*Y(J) - ...
                N(R(U))*v)./(N(R(i))+N(R(j))+N(R(U)));
    end
    J = [J i*(m-(i+1)/2)-m+j];
    Y(J) = []; % no need for the cluster information about j.

    % update m, N, R
    m = m-1;
    N(n+s) = N(R(i)) + N(R(j));
    R(i) = n+s;
    R(j:(n-1))=R((j+1):n);
end

%--------------------------------------------------------------------------

function changes = computeChanges(ind, adjprior_cq, adjprior_sp, ...
    indCqCounts, indSpCounts)
% Computes changes in log-marginal likelihood if individual ind is
% moved to another population
%
% Input:
% ind - the individual to be moved
% adjprior_cq & _sp  - adjpriors for cliques and separators
% indCqCounts, indSpCounts - counts for individual ind
%
% Output:
% changes - table of size 1*npops. changes(i) = difference in logml if
% ind is move to population i.

global CQ_COUNTS;   global CQ_SUMCOUNTS;
global SP_COUNTS;   global SP_SUMCOUNTS;
global PARTITION;   global POP_LOGML;
global LOGDIFF;

%npops = size(CQ_COUNTS,3); deleted by Lu Cheng, 24.20.2010
changes = LOGDIFF(ind,:);

i1 = PARTITION(ind);
i1_logml = POP_LOGML(i1);
changes(i1) = 0;

sumCq = uint16(sum(indCqCounts,1)); %n_loci * 1 vector
sumSp = uint16(sum(indSpCounts,1));

CQ_COUNTS(:,:,i1) = CQ_COUNTS(:,:,i1)-indCqCounts;
CQ_SUMCOUNTS(i1,:) = CQ_SUMCOUNTS(i1,:)-sumCq;
SP_COUNTS(:,:,i1) = SP_COUNTS(:,:,i1)-indSpCounts;
SP_SUMCOUNTS(i1,:) = SP_SUMCOUNTS(i1,:)-sumSp;

new_i1_logml = computePopulationLogml(i1, adjprior_cq, adjprior_sp);

CQ_COUNTS(:,:,i1) = CQ_COUNTS(:,:,i1)+indCqCounts;
CQ_SUMCOUNTS(i1,:) = CQ_SUMCOUNTS(i1,:)+sumCq;
SP_COUNTS(:,:,i1) = SP_COUNTS(:,:,i1)+indSpCounts;
SP_SUMCOUNTS(i1,:) = SP_SUMCOUNTS(i1,:)+sumSp;

i2 = find(changes==-Inf);
% i2 = setdiff(i2,i1);  % deleted by Lu Cheng, 24.02.2010, since changes(i1) = 0
i2_logml = POP_LOGML(i2);

ni2 = length(i2);

CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2)+repmat(indCqCounts, [1 1 ni2]);
CQ_SUMCOUNTS(i2,:) = CQ_SUMCOUNTS(i2,:)+repmat(sumCq,[ni2 1]);
SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2)+repmat(indSpCounts, [1 1 ni2]);
SP_SUMCOUNTS(i2,:) = SP_SUMCOUNTS(i2,:) + repmat(sumSp,[ni2 1]);

new_i2_logml = computePopulationLogml(i2, adjprior_cq, adjprior_sp);

CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2)-repmat(indCqCounts, [1 1 ni2]);
CQ_SUMCOUNTS(i2,:) = CQ_SUMCOUNTS(i2,:)-repmat(sumCq,[ni2 1]);
SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2)-repmat(indSpCounts, [1 1 ni2]);
SP_SUMCOUNTS(i2,:) = SP_SUMCOUNTS(i2,:) - repmat(sumSp,[ni2 1]);

changes(i2) = new_i1_logml - i1_logml ...
    + new_i2_logml - i2_logml;
LOGDIFF(ind,:) = changes;

%------------------------------------------------------------------------------------

function changes = computeChanges2(i1, adjprior_cq, adjprior_sp)
% Computes changes in log marginal likelihood if population i1 is combined
% with another population
%
% Input:
% i1 - the population to be combined
% adjprior_cq & _sp  - adjpriors for cliques and separators
%
% Output:
% changes - table of size 1*npops. changes(i) = difference in logml if
% i1 is combined with population i.

global CQ_COUNTS;   global CQ_SUMCOUNTS;
global SP_COUNTS;   global SP_SUMCOUNTS;
global POP_LOGML;

npops = size(CQ_COUNTS,3);
changes = zeros(npops,1);

i1_logml = POP_LOGML(i1);
indCqCounts = CQ_COUNTS(:,:,i1);
indSpCounts = SP_COUNTS(:,:,i1);
sumCq = uint16(sum(indCqCounts,1));
sumSp = uint16(sum(indSpCounts,1));

new_i1_logml = 0;

i2 = [1:i1-1 , i1+1:npops];
i2_logml = POP_LOGML(i2);

CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2)+repmat(indCqCounts, [1 1 npops-1]);
CQ_SUMCOUNTS(i2,:) = CQ_SUMCOUNTS(i2,:)+repmat(sumCq,[npops-1 1]);
SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2)+repmat(indSpCounts, [1 1 npops-1]);
SP_SUMCOUNTS(i2,:) = SP_SUMCOUNTS(i2,:)+ repmat(sumSp,[npops-1 1]);
% a = repmat(sumSp,[npops-1 1]);
% if ~any(sumSp)
%     a(:,[1:size(a,2)])=[];
% end
% SP_SUMCOUNTS(i2,:) = SP_SUMCOUNTS(i2,:)+ a ;


new_i2_logml = computePopulationLogml(i2, adjprior_cq, adjprior_sp);

CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2)-repmat(indCqCounts, [1 1 npops-1]);
CQ_SUMCOUNTS(i2,:) = CQ_SUMCOUNTS(i2,:)-repmat(sumCq,[npops-1 1]);
SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2)-repmat(indSpCounts, [1 1 npops-1]);
SP_SUMCOUNTS(i2,:) = SP_SUMCOUNTS(i2,:)- repmat(sumSp,[npops-1 1]);

changes(i2) = new_i1_logml - i1_logml + new_i2_logml - i2_logml;

%------------------------------------------------------------------------------------


function changes = computeChanges3(T2, inds2, i1, counts_cq, counts_sp, ...
    adjprior_cq, adjprior_sp)
% Computes changes in log marginal likelihood if subpopulation of i2 is
% moved to another population
%
% Input:
% T2 - partition of inds2 to subpopulations
% inds2 - individuals in population i1
% i2
% counts_cq, counts_sp - counts for individuals
%
% Output:
% changes - table of size length(unique(T2))*npops.
% changes(i,j) = difference in logml if subpopulation inds2(find(T2==i)) of
% i2 is moved to population j

global CQ_COUNTS;   global CQ_SUMCOUNTS;
global SP_COUNTS;   global SP_SUMCOUNTS;
global POP_LOGML;

npops = size(CQ_COUNTS,3);
npops2 = length(unique(T2));
changes = zeros(npops2,npops);

%cq_counts = CQ_COUNTS;
%sp_counts = SP_COUNTS;
%cq_sumcounts = CQ_SUMCOUNTS;
%sp_sumcounts = SP_SUMCOUNTS;


i1_logml = POP_LOGML(i1);

for pop2 = 1:npops2
    % inds = inds2(find(T2==pop2));
    inds = inds2(T2==pop2);
    ninds = length(inds);
    if ninds>0
        indCqCounts = uint16(sum(counts_cq(:,:,inds),3));
        indSpCounts = uint16(sum(counts_sp(:,:,inds),3));
        sumCq = uint16(sum(indCqCounts,1));
        sumSp = uint16(sum(indSpCounts,1));

        CQ_COUNTS(:,:,i1) = CQ_COUNTS(:,:,i1)-indCqCounts;
        CQ_SUMCOUNTS(i1,:) = CQ_SUMCOUNTS(i1,:)-sumCq;
        SP_COUNTS(:,:,i1) = SP_COUNTS(:,:,i1)-indSpCounts;
        SP_SUMCOUNTS(i1,:) = SP_SUMCOUNTS(i1,:)-sumSp;

        new_i1_logml = computePopulationLogml(i1, adjprior_cq, adjprior_sp);

        CQ_COUNTS(:,:,i1) = CQ_COUNTS(:,:,i1)+indCqCounts;
        CQ_SUMCOUNTS(i1,:) = CQ_SUMCOUNTS(i1,:)+sumCq;
        SP_COUNTS(:,:,i1) = SP_COUNTS(:,:,i1)+indSpCounts;
        SP_SUMCOUNTS(i1,:) = SP_SUMCOUNTS(i1,:)+sumSp;

        i2 = [1:i1-1 , i1+1:npops];
        i2_logml = POP_LOGML(i2)';

        CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2)+repmat(indCqCounts, [1 1 npops-1]);
        CQ_SUMCOUNTS(i2,:) = CQ_SUMCOUNTS(i2,:)+repmat(sumCq,[npops-1 1]);
        SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2)+repmat(indSpCounts, [1 1 npops-1]);
        SP_SUMCOUNTS(i2,:) = SP_SUMCOUNTS(i2,:)+ repmat(sumSp,[npops-1 1]);

        new_i2_logml = computePopulationLogml(i2, adjprior_cq, adjprior_sp)';

        CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2)-repmat(indCqCounts, [1 1 npops-1]);
        CQ_SUMCOUNTS(i2,:) = CQ_SUMCOUNTS(i2,:)-repmat(sumCq,[npops-1 1]);
        SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2)-repmat(indSpCounts, [1 1 npops-1]);
        SP_SUMCOUNTS(i2,:) = SP_SUMCOUNTS(i2,:)- repmat(sumSp,[npops-1 1]);

        changes(pop2,i2) = new_i1_logml - i1_logml + new_i2_logml - i2_logml;
    end
end

%--------------------------------------------------------------------------

function changes = computeChanges5(inds, i1, i2, counts_cq, counts_sp, ...
    adjprior_cq, adjprior_sp)
% Computes change in logml if individual of inds is moved between
% populations i1 and i2

global CQ_COUNTS;   global CQ_SUMCOUNTS;
global SP_COUNTS;   global SP_SUMCOUNTS;
global POP_LOGML;   global PARTITION;

ninds = length(inds);
changes = zeros(ninds,1);

i1_logml = POP_LOGML(i1);
i2_logml = POP_LOGML(i2);

for i = 1:ninds
    ind = inds(i);
    if PARTITION(ind)==i1
        pop1 = i1;  %from
        pop2 = i2;  %to
    else
        pop1 = i2;
        pop2 = i1;
    end
    indCqCounts = uint16(counts_cq(:,:,ind));
    indSpCounts = uint16(counts_sp(:,:,ind));
    sumCq = uint16(sum(indCqCounts,1));
    sumSp = uint16(sum(indSpCounts,1));

    CQ_COUNTS(:,:,pop1) = CQ_COUNTS(:,:,pop1)-indCqCounts;
    CQ_SUMCOUNTS(pop1,:) = CQ_SUMCOUNTS(pop1,:)-sumCq;
    SP_COUNTS(:,:,pop1) = SP_COUNTS(:,:,pop1)-indSpCounts;
    SP_SUMCOUNTS(pop1,:) = SP_SUMCOUNTS(pop1,:) - sumSp;

    CQ_COUNTS(:,:,pop2) = CQ_COUNTS(:,:,pop2)+indCqCounts;
    CQ_SUMCOUNTS(pop2,:) = CQ_SUMCOUNTS(pop2,:)+sumCq;
    SP_COUNTS(:,:,pop2) = SP_COUNTS(:,:,pop2)+indSpCounts;
    SP_SUMCOUNTS(pop2,:) = SP_SUMCOUNTS(pop2,:) + sumSp;

    new_logmls = computePopulationLogml([i1 i2], adjprior_cq, adjprior_sp);
    changes(i) = sum(new_logmls);

    CQ_COUNTS(:,:,pop1) = CQ_COUNTS(:,:,pop1)+indCqCounts;
    CQ_SUMCOUNTS(pop1,:) = CQ_SUMCOUNTS(pop1,:)+sumCq;
    SP_COUNTS(:,:,pop1) = SP_COUNTS(:,:,pop1)+indSpCounts;
    SP_SUMCOUNTS(pop1,:) = SP_SUMCOUNTS(pop1,:)+sumSp;
    CQ_COUNTS(:,:,pop2) = CQ_COUNTS(:,:,pop2)-indCqCounts;
    CQ_SUMCOUNTS(pop2,:) = CQ_SUMCOUNTS(pop2,:)-sumCq;
    SP_COUNTS(:,:,pop2) = SP_COUNTS(:,:,pop2)-indSpCounts;
    SP_SUMCOUNTS(pop2,:) = SP_SUMCOUNTS(pop2,:)-sumSp;
end

changes = changes - i1_logml - i2_logml;


%-------------------------------------------------------------------------------------


function updateGlobalVariables(ind, i2, indCqCounts, indSpCounts, ...
    adjprior_cq, adjprior_sp)
% Updates global variables when individual ind is moved to population i2

global CQ_COUNTS;   global CQ_SUMCOUNTS;
global SP_COUNTS;   global SP_SUMCOUNTS;
global PARTITION;   global POP_LOGML;
global LOGDIFF;

i1 = PARTITION(ind);
PARTITION(ind)=i2;

sumCq = uint16(sum(indCqCounts,1));
sumSp = uint16(sum(indSpCounts,1));

CQ_COUNTS(:,:,i1) = CQ_COUNTS(:,:,i1)-indCqCounts;
CQ_SUMCOUNTS(i1,:) = CQ_SUMCOUNTS(i1,:)-sumCq;
SP_COUNTS(:,:,i1) = SP_COUNTS(:,:,i1)-indSpCounts;
SP_SUMCOUNTS(i1,:) = SP_SUMCOUNTS(i1,:)-sumSp;

CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2)+indCqCounts;
CQ_SUMCOUNTS(i2,:) = CQ_SUMCOUNTS(i2,:)+sumCq;
SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2)+indSpCounts;
SP_SUMCOUNTS(i2,:) = SP_SUMCOUNTS(i2,:)+sumSp;

POP_LOGML([i1 i2]) = computePopulationLogml([i1 i2], adjprior_cq, adjprior_sp);

LOGDIFF(:,[i1 i2]) = -Inf;
%inx = [find(PARTITION==i1); find(PARTITION==i2)];
inx = (PARTITION==i1 | PARTITION==i2); % modified by Lu Cheng, 24.02.2010
LOGDIFF(inx,:) = -Inf;


%---------------------------------------------------------------------------------


function updateGlobalVariables2(i1, i2, adjprior_cq, adjprior_sp)
% Updates global variables when all individuals from population i1 are moved
% to population i2

global CQ_COUNTS;   global CQ_SUMCOUNTS;
global SP_COUNTS;   global SP_SUMCOUNTS;
global PARTITION;   global POP_LOGML;
global LOGDIFF;

% inds = find(PARTITION==i1);
% PARTITION(inds) = i2;
PARTITION(PARTITION==i1) = i2;

CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2)+CQ_COUNTS(:,:,i1);
CQ_SUMCOUNTS(i2,:) = CQ_SUMCOUNTS(i2,:)+CQ_SUMCOUNTS(i1,:);
SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2)+SP_COUNTS(:,:,i1);
SP_SUMCOUNTS(i2,:) = SP_SUMCOUNTS(i2,:)+SP_SUMCOUNTS(i1,:);

CQ_COUNTS(:,:,i1) = 0;
CQ_SUMCOUNTS(i1,:) = 0;
SP_COUNTS(:,:,i1) = 0;
SP_SUMCOUNTS(i1,:) = 0;

POP_LOGML(i1) = 0;
POP_LOGML(i2) = computePopulationLogml(i2, adjprior_cq, adjprior_sp);

LOGDIFF(:,[i1 i2]) = -Inf;
%inx = [find(PARTITION==i1); find(PARTITION==i2)];
inx = (PARTITION==i1 | PARTITION==i2); % modified by Lu Cheng, 24.02.2010
LOGDIFF(inx,:) = -Inf;

%------------------------------------------------------------------------------------


function updateGlobalVariables3(muuttuvat, i2, indCqCounts, indSpCounts, ...
    adjprior_cq, adjprior_sp)
% Updates global variables when individuals muuttuvat are moved to
% population i2

global CQ_COUNTS;   global CQ_SUMCOUNTS;
global SP_COUNTS;   global SP_SUMCOUNTS;
global PARTITION;   global POP_LOGML;
global LOGDIFF;

i1 = PARTITION(muuttuvat(1));   %Why only one individual is moved? Lu Cheng, 24.02.2010
PARTITION(muuttuvat) = i2;

sumCq = uint16(sum(indCqCounts,1));
sumSp = uint16(sum(indSpCounts,1));

CQ_COUNTS(:,:,i1) = CQ_COUNTS(:,:,i1)-indCqCounts;
CQ_SUMCOUNTS(i1,:) = CQ_SUMCOUNTS(i1,:)-sumCq;
SP_COUNTS(:,:,i1) = SP_COUNTS(:,:,i1)-indSpCounts;
SP_SUMCOUNTS(i1,:) = SP_SUMCOUNTS(i1,:)-sumSp;

CQ_COUNTS(:,:,i2) = CQ_COUNTS(:,:,i2)+indCqCounts;
CQ_SUMCOUNTS(i2,:) = CQ_SUMCOUNTS(i2,:)+sumCq;
SP_COUNTS(:,:,i2) = SP_COUNTS(:,:,i2)+indSpCounts;
SP_SUMCOUNTS(i2,:) = SP_SUMCOUNTS(i2,:)+sumSp;

POP_LOGML([i1 i2]) = computePopulationLogml([i1 i2], adjprior_cq, adjprior_sp);

LOGDIFF(:,[i1 i2]) = -Inf;
inx = [find(PARTITION==i1); find(PARTITION==i2)];
LOGDIFF(inx,:) = -Inf;

%----------------------------------------------------------------------


function inds = returnInOrder(inds, pop, counts_cq, counts_sp, ...
    adjprior_cq, adjprior_sp)
% Returns individuals inds in order according to the change in the logml if
% they are moved out of the population pop

global CQ_COUNTS;   global CQ_SUMCOUNTS;
global SP_COUNTS;   global SP_SUMCOUNTS;

ninds = length(inds);
apuTaulu = [inds, zeros(ninds,1)];

for i=1:ninds
    ind = inds(i);
    indCqCounts = uint16(counts_cq(:,:,ind));
    indSpCounts = uint16(counts_sp(:,:,ind));
    sumCq = uint16(sum(indCqCounts,1));
    sumSp = uint16(sum(indSpCounts,1));

    CQ_COUNTS(:,:,pop) = CQ_COUNTS(:,:,pop)-indCqCounts;
    CQ_SUMCOUNTS(pop,:) = CQ_SUMCOUNTS(pop,:)-sumCq;
    SP_COUNTS(:,:,pop) = SP_COUNTS(:,:,pop)-indSpCounts;
    SP_SUMCOUNTS(pop,:) = SP_SUMCOUNTS(pop,:)-sumSp;

    apuTaulu(i, 2) = computePopulationLogml(pop, adjprior_cq, adjprior_sp);

    CQ_COUNTS(:,:,pop) = CQ_COUNTS(:,:,pop)+indCqCounts;
    CQ_SUMCOUNTS(pop,:) = CQ_SUMCOUNTS(pop,:)+sumCq;
    SP_COUNTS(:,:,pop) = SP_COUNTS(:,:,pop)+indSpCounts;
    SP_SUMCOUNTS(pop,:) = SP_SUMCOUNTS(pop,:)+sumSp;
end
apuTaulu = sortrows(apuTaulu,2);
inds = apuTaulu(ninds:-1:1,1);


%--------------------------------------------------------------------------

function clearGlobalVars

global CQ_COUNTS; CQ_COUNTS = [];
global CQ_SUMCOUNTS; CQ_SUMCOUNTS = [];
global SP_COUNTS; SP_COUNTS = [];
global SP_SUMCOUNTS; SP_SUMCOUNTS = [];
global PARTITION; PARTITION = [];
global POP_LOGML; POP_LOGML = [];
global LOGDIFF; LOGDIFF = [];

global TRAIN_CQ_COUNTS;    TRAIN_CQ_COUNTS=[];
global TRAIN_CQ_SUMCOUNTS; TRAIN_CQ_SUMCOUNTS = [];
global TRAIN_SP_COUNTS;    TRAIN_SP_COUNTS = [];
global TRAIN_SP_SUMCOUNTS; TRAIN_SP_SUMCOUNTS = [];


%--------------------------------------------------------------------------

function new_npops = removeEmptyPops(n_train_clusters, npops)
% Removes empty pops from all global COUNTS variables.
% Updates PARTITION and npops

% This function is modified by Lu Cheng, 08.03.2010
% n_train_clusters: number of clusters in the training data
% npops: the user defined max potential population number

global CQ_COUNTS;
global CQ_SUMCOUNTS;
global SP_COUNTS;
global SP_SUMCOUNTS;
global PARTITION;
global LOGDIFF;

global TRAIN_CQ_COUNTS;
global TRAIN_CQ_SUMCOUNTS;
global TRAIN_SP_COUNTS;
global TRAIN_SP_SUMCOUNTS;

%---------ORIGINAL CODES START------------------%
% notEmpty = find(any(CQ_SUMCOUNTS,2));
% CQ_COUNTS = CQ_COUNTS(:,:,notEmpty);
% CQ_SUMCOUNTS = CQ_SUMCOUNTS(notEmpty,:);
% SP_COUNTS = SP_COUNTS(:,:,notEmpty);
% SP_SUMCOUNTS = SP_SUMCOUNTS(notEmpty,:);
% LOGDIFF = LOGDIFF(:,notEmpty);
%
% for n=1:length(notEmpty)
%     PARTITION(PARTITION==notEmpty(n)) = n;
% end
% npops = length(notEmpty);
%---------ORIGINAL CODES END---------------------%

pops = unique(PARTITION);
pops = pops(pops>n_train_clusters);
if isempty(pops) && npops>n_train_clusters
    CQ_COUNTS(:,:,n_train_clusters+1:npops) = [];
    CQ_SUMCOUNTS(n_train_clusters+1:npops,:) = [];
    SP_COUNTS(:,:,n_train_clusters+1:npops) = [];
    SP_SUMCOUNTS(n_train_clusters+1:npops,:) = [];

    TRAIN_CQ_COUNTS(:,:,n_train_clusters+1:npops) = [];
    TRAIN_CQ_SUMCOUNTS(n_train_clusters+1:npops,:) = [];
    TRAIN_SP_COUNTS(:,:,n_train_clusters+1:npops) = [];
    TRAIN_SP_SUMCOUNTS(n_train_clusters+1:npops,:) = [];

    LOGDIFF(:,n_train_clusters+1:npops) = [];
    new_npops = n_train_clusters;
    return
elseif isempty(pops)
    % do nothing

    new_npops = npops;
    return
else
    n_nonempty = length(pops); %the pops should be in ascending order
    for i=1:n_nonempty
        PARTITION(PARTITION==pops(i))= n_train_clusters+i;
        CQ_COUNTS(:,:,n_train_clusters+i) = CQ_COUNTS(:,:,pops(i));
        CQ_SUMCOUNTS(n_train_clusters+i,:) = CQ_SUMCOUNTS(pops(i),:);
        SP_COUNTS(:,:,n_train_clusters+i) = SP_COUNTS(:,:,pops(i));
        SP_SUMCOUNTS(n_train_clusters+i,:) = SP_SUMCOUNTS(pops(i),:);

        LOGDIFF(:,n_train_clusters+i) = LOGDIFF(:,pops(i));
    end

    CQ_COUNTS(:,:,n_train_clusters+n_nonempty+1:npops) = [];
    CQ_SUMCOUNTS(n_train_clusters+n_nonempty+1:npops,:) = [];
    SP_COUNTS(:,:,n_train_clusters+n_nonempty+1:npops) = [];
    SP_SUMCOUNTS(n_train_clusters+n_nonempty+1:npops,:) = [];

    TRAIN_CQ_COUNTS(:,:,n_train_clusters+n_nonempty+1:npops) = [];
    TRAIN_CQ_SUMCOUNTS(n_train_clusters+n_nonempty+1:npops,:) = [];
    TRAIN_SP_COUNTS(:,:,n_train_clusters+n_nonempty+1:npops) = [];
    TRAIN_SP_SUMCOUNTS(n_train_clusters+n_nonempty+1:npops,:) = [];

    LOGDIFF(:,n_train_clusters+n_nonempty+1:npops) = [];

    new_npops = n_train_clusters+n_nonempty;
end



%--------------------------------------------------------------------------

function [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex, n_train_pops)
% Tiedet��n, ett� annettu logml on isompi kuin huonoin arvo
% partitionSummary taulukossa. Jos partitionSummary:ss� ei viel� ole
% annettua logml arvoa, niin lis�t��n worstIndex:in kohtaan uusi logml ja
% nykyist� partitiota vastaava nclusters:in arvo. Muutoin ei tehd� mit��n.

% global PARTITION;
% apu = isempty(find(abs(partitionSummary(:,2)-logml)<1e-5,1));
% if apu
%     % Nyt l�ydetty partitio ei ole viel� kirjattuna summaryyn.
%     npops = length(unique(PARTITION));
%     partitionSummary(worstIndex,1) = npops;
%     partitionSummary(worstIndex,2) = logml;
%     added = 1;
% else
%     added = 0;
% end

% Modified by Lu Cheng, 11.03.2010
% The population number should be greater than or equal to the training
% population number
global PARTITION;
apu = isempty(find(abs(partitionSummary(:,2)-logml)<1e-5,1));
if apu
    npops = sum(unique(PARTITION)>n_train_pops)+n_train_pops;
    partitionSummary(worstIndex,1) = npops;
    partitionSummary(worstIndex,2) = logml;
    added = 1;
else
    added = 0;
end

%--------------------------------------------------------------------------

function [counts, sumcounts] = initialCounts(ind_counts)

global PARTITION;

pops = unique(PARTITION);
npops = max(pops);

counts = zeros(size(ind_counts,1), size(ind_counts,2), npops,'uint16');
sumcounts = zeros(npops, size(ind_counts,2),'uint16');

for i = 1:npops
    inds = (PARTITION == i);
    counts(:,:,i) = sum(ind_counts(:,:,inds), 3);
    sumcounts(i,:) = sum(counts(:,:,i),1);
end

%--------------------------------------------------------------------------

function [counts, sumcounts] = calCounts(ind_counts, npops)
% This function calculates the counts for each cluster from the sample data
% 'ind_counts', which are the counts for each individual
% The input 'npops' is a user defined value, which means the max number of
% potential populations in the data. Usually it is bigger than the number
% of populations given by 'PARTITION'.

% modified by Lu Cheng, 08.03.2010

global PARTITION;  % n_samples * 1 vector, indicate the cluster_id for each sample

pops = unique(PARTITION);
pops_len = length(pops);

[n_alle n_loci n_inds] = size(ind_counts);
counts = zeros(n_alle, n_loci, npops,'uint16');
sumcounts = zeros(npops, n_loci,'uint16');

for i = 1:pops_len
    cluster = pops(i);
    inds = (PARTITION == cluster);
    counts(:,:,cluster) = sum(ind_counts(:,:,inds), 3);
    sumcounts(cluster,:) = squeeze(sum(counts(:,:,cluster),1));
end

%--------------------------------------------------------------------------

function logml = computeLogml(adjprior_cq, adjprior_sp)

% calculate the sum logml for all the populations

global CQ_COUNTS;   global CQ_SUMCOUNTS;
global SP_COUNTS;   global SP_SUMCOUNTS;

global TRAIN_CQ_COUNTS;
global TRAIN_SP_COUNTS;
global TRAIN_CQ_SUMCOUNTS;
global TRAIN_SP_SUMCOUNTS;

npops = size(TRAIN_CQ_COUNTS,3);                     % Modified by Lu Cheng, 08.03.2010
cq_counts = TRAIN_CQ_COUNTS + double(CQ_COUNTS);
cq_sumcounts = TRAIN_CQ_SUMCOUNTS + double(CQ_SUMCOUNTS);
sp_counts = TRAIN_SP_COUNTS + double(SP_COUNTS);
sp_sumcounts = TRAIN_SP_SUMCOUNTS + double(SP_SUMCOUNTS);

cq_logml = sum(sum(sum(gammaln(cq_counts+repmat(adjprior_cq,[1 1 npops]))))) ...
    - npops*sum(sum(gammaln(adjprior_cq))) - ...
    sum(sum(gammaln(1+cq_sumcounts)));

sp_logml = sum(sum(sum(gammaln(sp_counts+repmat(adjprior_sp,[1 1 npops]))))) ...
    - npops*sum(sum(gammaln(adjprior_sp))) - ...
    sum(sum(gammaln(1+sp_sumcounts)));

logml = cq_logml - sp_logml;
clear cq_counts cq_sumcounts sp_counts sp_sumcounts;

%--------------------------------------------------------------------------

function popLogml = computePopulationLogml(pops, adjprior_cq, adjprior_sp)
% Palauttaa length(pops)*1 taulukon, jossa on laskettu korikohtaiset
% logml:t koreille, jotka on m��ritelty pops-muuttujalla.

% calculate logml for each population given by pops
% returns a length(pops)*1 vector

global CQ_COUNTS;   global CQ_SUMCOUNTS;
global SP_COUNTS;   global SP_SUMCOUNTS;

global TRAIN_CQ_COUNTS;
global TRAIN_SP_COUNTS;
global TRAIN_CQ_SUMCOUNTS;
global TRAIN_SP_SUMCOUNTS;

%npops = size(TRAIN_CQ_COUNTS,3);                     % Modified by Lu Cheng, 08.03.2010
cq_counts = TRAIN_CQ_COUNTS + double(CQ_COUNTS);
cq_sumcounts = TRAIN_CQ_SUMCOUNTS + double(CQ_SUMCOUNTS);
sp_counts = TRAIN_SP_COUNTS + double(SP_COUNTS);
sp_sumcounts = TRAIN_SP_SUMCOUNTS + double(SP_SUMCOUNTS);

nall_cq = size(CQ_COUNTS,1);
nall_sp = size(SP_COUNTS,1);
ncliq = size(CQ_COUNTS,2);
nsep = size(SP_COUNTS, 2);

z = length(pops);

popLogml_cq = ...
    squeeze(sum(sum(reshape(...
    gammaln(repmat(adjprior_cq,[1 1 z]) + cq_counts(:,:,pops)) ...
    ,[nall_cq ncliq z]),1),2)) - sum(gammaln(1+cq_sumcounts(pops,:)),2) - ...
    sum(sum(gammaln(adjprior_cq)));

popLogml_sp = ...
    squeeze(sum(sum(reshape(...
    gammaln(repmat(adjprior_sp,[1 1 z]) + sp_counts(:,:,pops)) ...
    ,[nall_sp nsep z]),1),2)) - sum(gammaln(1+sp_sumcounts(pops,:)),2) - ...
    sum(sum(gammaln(adjprior_sp)));

popLogml = popLogml_cq - popLogml_sp;
clear cq_counts cq_sumcounts sp_counts sp_sumcounts;

%-------------------------------------------------------------------




%--------------------------------------------------------------
function newline = takeLine(description,width)
%Returns one line from the description: line ends to the first
%space after width:th mark.
% newLine = description(1:width);
n = width+1;
while ~isspace(description(n)) && n<length(description)
    n = n+1;
end;
newline = description(1:n);


function dispLine
disp('---------------------------------------------------');

function dispCancel
disp('** CANCELLED');

function num2 = omaRound(num)
% Py�rist�� luvun num 1 desimaalin tarkkuuteen
num = num*10;
num = round(num);
num2 = num/10;

%---------------------------------------------------------
function mjono = logml2String(logml)
% Palauttaa logml:n string-esityksen.

mjono = '       ';
if abs(logml)<10000
    %Ei tarvita e-muotoa
    mjono(7) = palautaYks(abs(logml),-1);
    mjono(6) = '.';
    mjono(5) = palautaYks(abs(logml),0);
    mjono(4) = palautaYks(abs(logml),1);
    mjono(3) = palautaYks(abs(logml),2);
    mjono(2) = palautaYks(abs(logml),3);
    pointer = 2;
    while mjono(pointer)=='0' && pointer<7
        mjono(pointer) = ' ';
        pointer=pointer+1;
    end
    if logml<0
        mjono(pointer-1) = '-';
    end
else
    suurinYks = 4;
    while abs(logml)/(10^(suurinYks+1)) >= 1
        suurinYks = suurinYks+1;
    end
    if suurinYks<10
        mjono(7) = num2str(suurinYks);
        mjono(6) = 'e';
        mjono(5) = palautaYks(abs(logml),suurinYks-1);
        mjono(4) = '.';
        mjono(3) = palautaYks(abs(logml),suurinYks);
        if logml<0
            mjono(2) = '-';
        end
    elseif suurinYks>=10
        mjono(6:7) = num2str(suurinYks);
        mjono(5) = 'e';
        mjono(4) = palautaYks(abs(logml),suurinYks-1);
        mjono(3) = '.';
        mjono(2) = palautaYks(abs(logml),suurinYks);
        if logml<0
            mjono(1) = '-';
        end
    end
end

function digit = palautaYks(num,yks)
% palauttaa luvun num 10^yks termin kertoimen
% string:in�
% yks t�ytyy olla kokonaisluku, joka on
% v�hint��n -1:n suuruinen. Pienemmill�
% luvuilla tapahtuu jokin py�ristysvirhe.

if yks>=0
    digit = rem(num, 10^(yks+1));
    digit = floor(digit/(10^yks));
else
    digit = num*10;
    digit = floor(rem(digit,10));
end
digit = num2str(digit);


%--------------------------------------------------------------------------

function [emptyPop, pops] = findEmptyPop(npops)
% Palauttaa ensimm�isen tyhj�n populaation indeksin. Jos tyhji�
% populaatioita ei ole, palauttaa -1:n.

global PARTITION;
pops = unique(PARTITION)';
if (length(pops) ==npops)
    emptyPop = -1;
else
    popDiff = diff([0 pops npops+1]);
    emptyPop = min(find(popDiff > 1));
end

%--------------------------------------------------------------------------
