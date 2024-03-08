function clearGlobalVars

global COUNTS; COUNTS = [];
global SUMCOUNTS; SUMCOUNTS = [];
global PARTITION; PARTITION = [];
global POP_LOGML; POP_LOGML = [];
global LOGDIFF; LOGDIFF = [];

%--------------------------------------------------------------------------


function Z = linkage(Y, method)
[k, n] = size(Y);
m = (1+sqrt(1+8*n))/2;
if k ~= 1 | m ~= fix(m)
  error('The first input has to match the output of the PDIST function in size.');
end
if nargin == 1 % set default switch to be 'co'
   method = 'co';
end
method = lower(method(1:2)); % simplify the switch string.
monotonic = 1;
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


%-----------------------------------------------------------------------


function [sumcounts, counts, logml] = ...
    initialPopCounts(data, npops, rows, noalle, adjprior)

nloci=size(data,2);
counts = zeros(max(noalle),nloci,npops);
sumcounts = zeros(npops,nloci);

for i=1:npops
    for j=1:nloci
        i_rivit = rows(i,1):rows(i,2);
        havainnotLokuksessa = find(data(i_rivit,j)>=0);
        sumcounts(i,j) = length(havainnotLokuksessa);
        for k=1:noalle(j)
            alleleCode = k;
            N_ijk = length(find(data(i_rivit,j)==alleleCode));
            counts(k,j,i) = N_ijk;
        end
    end
end

logml = laskeLoggis(counts, sumcounts, adjprior);

%-----------------------------------------------------------------------


function loggis = laskeLoggis(counts, sumcounts, adjprior)
npops = size(counts,3);

logml2 = sum(sum(sum(gammaln(counts+repmat(adjprior,[1 1 npops]))))) ...
    - npops*sum(sum(gammaln(adjprior))) - ...
    sum(sum(gammaln(1+sumcounts)));
loggis = logml2;


%------------------------------------------------------------------------------------

function [muutokset, diffInCounts] = ...
    laskeMuutokset(ind, globalRows, da20ta, adjprior, priorTerm)
% Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mikה olisi
% muutos logml:ssה, mikהli yksilצ ind siirretההn koriin i.
% diffInCounts on poistettava COUNTS:in siivusta i1 ja lisהttהvה
% COUNTS:in siivuun i2, mikהli muutos toteutetaan.
%
% Lisהys 25.9.2007:
% Otettu kהyttצצn globaali muuttuja LOGDIFF, johon on tallennettu muutokset
% logml:ssה siirrettהessה yksilצitה toisiin populaatioihin.

global COUNTS;      global SUMCOUNTS;
global PARTITION;   global POP_LOGML;
global LOGDIFF;

npops = size(COUNTS,3);
muutokset = LOGDIFF(ind,:);

i1 = PARTITION(ind);
i1_logml = POP_LOGML(i1);
muutokset(i1) = 0;

rows = globalRows(ind,1):globalRows(ind,2);
diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
diffInSumCounts = sum(diffInCounts);

COUNTS(:,:,i1) = COUNTS(:,:,i1)-diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)-diffInSumCounts;
new_i1_logml = computePopulationLogml(i1, adjprior, priorTerm);
COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;

i2 = find(muutokset==-Inf);     % Etsitההn populaatiot jotka muuttuneet viime kerran jהlkeen.
i2 = setdiff(i2,i1);
i2_logml = POP_LOGML(i2);

ni2 = length(i2);

COUNTS(:,:,i2) = COUNTS(:,:,i2)+repmat(diffInCounts, [1 1 ni2]);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)+repmat(diffInSumCounts,[ni2 1]);
new_i2_logml = computePopulationLogml(i2, adjprior, priorTerm);
COUNTS(:,:,i2) = COUNTS(:,:,i2)-repmat(diffInCounts, [1 1 ni2]);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)-repmat(diffInSumCounts,[ni2 1]);

muutokset(i2) = new_i1_logml - i1_logml ...
    + new_i2_logml - i2_logml;
LOGDIFF(ind,:) = muutokset;


%----------------------------------------------------------------------


function diffInCounts = computeDiffInCounts(rows, max_noalle, nloci, data)
% Muodostaa max_noalle*nloci taulukon, jossa on niiden alleelien
% lukumההrהt (vastaavasti kuin COUNTS:issa), jotka ovat data:n
% riveillה rows. rows pitהה olla vaakavektori.

diffInCounts = zeros(max_noalle, nloci);
for i=rows
    row = data(i,:);
    notEmpty = find(row>=0);

    if length(notEmpty)>0
        diffInCounts(row(notEmpty) + (notEmpty-1)*max_noalle) = ...
            diffInCounts(row(notEmpty) + (notEmpty-1)*max_noalle) + 1;
    end
end

%------------------------------------------------------------------------


%-------------------------------------------------------------------------------------


function updateGlobalVariables(ind, i2, diffInCounts, ...
    adjprior, priorTerm)
% Suorittaa globaalien muuttujien muutokset, kun yksilצ ind
% on siirretההn koriin i2.

global PARTITION;
global COUNTS;
global SUMCOUNTS;
global POP_LOGML;
global LOGDIFF;

i1 = PARTITION(ind);
PARTITION(ind)=i2;

COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - sum(diffInCounts);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + sum(diffInCounts);

POP_LOGML([i1 i2]) = computePopulationLogml([i1 i2], adjprior, priorTerm);

LOGDIFF(:,[i1 i2]) = -Inf;
inx = [find(PARTITION==i1); find(PARTITION==i2)];
LOGDIFF(inx,:) = -Inf;


%--------------------------------------------------------------------------
%--

%------------------------------------------------------------------------------------


function [muutokset, diffInCounts] = laskeMuutokset2( ...
    i1, globalRows, data, adjprior, priorTerm);
% Palauttaa npops*1 taulun, jossa i:s alkio kertoo, mikה olisi
% muutos logml:ssה, mikהli korin i1 kaikki yksilצt siirretההn
% koriin i.

global COUNTS;      global SUMCOUNTS;
global PARTITION;   global POP_LOGML;
npops = size(COUNTS,3);
muutokset = zeros(npops,1);

i1_logml = POP_LOGML(i1);

inds = find(PARTITION==i1);
ninds = length(inds);

if ninds==0
    diffInCounts = zeros(size(COUNTS,1), size(COUNTS,2));
    return;
end

rows = [];
for i = 1:ninds
    ind = inds(i);
    lisa = globalRows(ind,1):globalRows(ind,2);
    rows = [rows; lisa'];
    %rows = [rows; globalRows{ind}'];
end

diffInCounts = computeDiffInCounts(rows', size(COUNTS,1), size(COUNTS,2), data);
diffInSumCounts = sum(diffInCounts);

COUNTS(:,:,i1) = COUNTS(:,:,i1)-diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)-diffInSumCounts;
new_i1_logml = computePopulationLogml(i1, adjprior, priorTerm);
COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;

i2 = [1:i1-1 , i1+1:npops];
i2_logml = POP_LOGML(i2);

COUNTS(:,:,i2) = COUNTS(:,:,i2)+repmat(diffInCounts, [1 1 npops-1]);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)+repmat(diffInSumCounts,[npops-1 1]);
new_i2_logml = computePopulationLogml(i2, adjprior, priorTerm);
COUNTS(:,:,i2) = COUNTS(:,:,i2)-repmat(diffInCounts, [1 1 npops-1]);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)-repmat(diffInSumCounts,[npops-1 1]);

muutokset(i2) = new_i1_logml - i1_logml ...
    + new_i2_logml - i2_logml;


%---------------------------------------------------------------------------------


function updateGlobalVariables2( ...
    i1, i2, diffInCounts, adjprior, priorTerm);
% Suorittaa globaalien muuttujien muutokset, kun kaikki
% korissa i1 olevat yksilצt siirretההn koriin i2.

global PARTITION;
global COUNTS;
global SUMCOUNTS;
global POP_LOGML;
global LOGDIFF;

inds = find(PARTITION==i1);
PARTITION(inds) = i2;

COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - sum(diffInCounts);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + sum(diffInCounts);

POP_LOGML(i1) = 0;
POP_LOGML(i2) = computePopulationLogml(i2, adjprior, priorTerm);

LOGDIFF(:,[i1 i2]) = -Inf;
inx = [find(PARTITION==i1); find(PARTITION==i2)];
LOGDIFF(inx,:) = -Inf;


%--------------------------------------------------------------------------
%----

function muutokset = laskeMuutokset3(T2, inds2, globalRows, ...
    data, adjprior, priorTerm, i1)
% Palauttaa length(unique(T2))*npops taulun, jossa (i,j):s alkio
% kertoo, mikה olisi muutos logml:ssה, jos populaation i1 osapopulaatio
% inds2(find(T2==i)) siirretההn koriin j.

global COUNTS;      global SUMCOUNTS;
global PARTITION;   global POP_LOGML;
npops = size(COUNTS,3);
npops2 = length(unique(T2));
muutokset = zeros(npops2, npops);

i1_logml = POP_LOGML(i1);
for pop2 = 1:npops2
    inds = inds2(find(T2==pop2));
    ninds = length(inds);
    if ninds>0
        rows = [];
        for i = 1:ninds
            ind = inds(i);
            lisa = globalRows(ind,1):globalRows(ind,2);
            rows = [rows; lisa'];
            %rows = [rows; globalRows{ind}'];
        end
        diffInCounts = computeDiffInCounts(rows', size(COUNTS,1), size(COUNTS,2), data);
        diffInSumCounts = sum(diffInCounts);

        COUNTS(:,:,i1) = COUNTS(:,:,i1)-diffInCounts;
        SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)-diffInSumCounts;
        new_i1_logml = computePopulationLogml(i1, adjprior, priorTerm);
        COUNTS(:,:,i1) = COUNTS(:,:,i1)+diffInCounts;
        SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:)+diffInSumCounts;

        i2 = [1:i1-1 , i1+1:npops];
        i2_logml = POP_LOGML(i2)';

        COUNTS(:,:,i2) = COUNTS(:,:,i2)+repmat(diffInCounts, [1 1 npops-1]);
        SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)+repmat(diffInSumCounts,[npops-1 1]);
        new_i2_logml = computePopulationLogml(i2, adjprior, priorTerm)';
        COUNTS(:,:,i2) = COUNTS(:,:,i2)-repmat(diffInCounts, [1 1 npops-1]);
        SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:)-repmat(diffInSumCounts,[npops-1 1]);

        muutokset(pop2,i2) = new_i1_logml - i1_logml ...
            + new_i2_logml - i2_logml;
    end
end

%------------------------------------------------------------------------------------

function muutokset = laskeMuutokset5(inds, globalRows, data, adjprior, ...
    priorTerm, i1, i2)

% Palauttaa length(inds)*1 taulun, jossa i:s alkio kertoo, mikה olisi
% muutos logml:ssה, mikהli yksilצ i vaihtaisi koria i1:n ja i2:n vהlillה.

global COUNTS;      global SUMCOUNTS;
global PARTITION;   global POP_LOGML;

ninds = length(inds);
muutokset = zeros(ninds,1);

i1_logml = POP_LOGML(i1);
i2_logml = POP_LOGML(i2);

for i = 1:ninds
    ind = inds(i);
    if PARTITION(ind)==i1
        pop1 = i1;  %mistה
        pop2 = i2;  %mihin
    else
        pop1 = i2;
        pop2 = i1;
    end
    rows = globalRows(ind,1):globalRows(ind,2);
    diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
    diffInSumCounts = sum(diffInCounts);

    COUNTS(:,:,pop1) = COUNTS(:,:,pop1)-diffInCounts;
    SUMCOUNTS(pop1,:) = SUMCOUNTS(pop1,:)-diffInSumCounts;
    COUNTS(:,:,pop2) = COUNTS(:,:,pop2)+diffInCounts;
    SUMCOUNTS(pop2,:) = SUMCOUNTS(pop2,:)+diffInSumCounts;

    new_logmls = computePopulationLogml([i1 i2], adjprior, priorTerm);
    muutokset(i) = sum(new_logmls);

    COUNTS(:,:,pop1) = COUNTS(:,:,pop1)+diffInCounts;
    SUMCOUNTS(pop1,:) = SUMCOUNTS(pop1,:)+diffInSumCounts;
    COUNTS(:,:,pop2) = COUNTS(:,:,pop2)-diffInCounts;
    SUMCOUNTS(pop2,:) = SUMCOUNTS(pop2,:)-diffInSumCounts;
end

muutokset = muutokset - i1_logml - i2_logml;

%------------------------------------------------------------------------------------


function updateGlobalVariables3(muuttuvat, diffInCounts, ...
    adjprior, priorTerm, i2);
% Suorittaa globaalien muuttujien pהivitykset, kun yksilצt 'muuttuvat'
% siirretההn koriin i2. Ennen siirtoa yksilצiden on kuuluttava samaan
% koriin.

global PARTITION;
global COUNTS;
global SUMCOUNTS;
global POP_LOGML;
global LOGDIFF;

i1 = PARTITION(muuttuvat(1));
PARTITION(muuttuvat) = i2;

COUNTS(:,:,i1) = COUNTS(:,:,i1) - diffInCounts;
COUNTS(:,:,i2) = COUNTS(:,:,i2) + diffInCounts;
SUMCOUNTS(i1,:) = SUMCOUNTS(i1,:) - sum(diffInCounts);
SUMCOUNTS(i2,:) = SUMCOUNTS(i2,:) + sum(diffInCounts);

POP_LOGML([i1 i2]) = computePopulationLogml([i1 i2], adjprior, priorTerm);

LOGDIFF(:,[i1 i2]) = -Inf;
inx = [find(PARTITION==i1); find(PARTITION==i2)];
LOGDIFF(inx,:) = -Inf;


%----------------------------------------------------------------------------


function dist2 = laskeOsaDist(inds2, dist, ninds)
% Muodostaa dist vektorista osavektorin, joka sisהltהה yksilצiden inds2
% vהliset etהisyydet. ninds=kaikkien yksilצiden lukumההrה.

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


%-----------------------------------------------------------------------------------


function npops = poistaTyhjatPopulaatiot(npops)
% Poistaa tyhjentyneet populaatiot COUNTS:ista ja
% SUMCOUNTS:ista. Pהivittהה npops:in ja PARTITION:in.

global COUNTS;
global SUMCOUNTS;
global PARTITION;
global LOGDIFF;

notEmpty = find(any(SUMCOUNTS,2));
COUNTS = COUNTS(:,:,notEmpty);
SUMCOUNTS = SUMCOUNTS(notEmpty,:);
LOGDIFF = LOGDIFF(:,notEmpty);

for n=1:length(notEmpty)
    apu = find(PARTITION==notEmpty(n));
    PARTITION(apu)=n;
end
npops = length(notEmpty);


%---------------------------------------------------------------


function dispLine;
disp('---------------------------------------------------');

%--------------------------------------------------------------

function num2 = omaRound(num)
% Pyצristהה luvun num 1 desimaalin tarkkuuteen
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
    while mjono(pointer)=='0' & pointer<7
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
% string:inה
% yks tהytyy olla kokonaisluku, joka on
% vהhintההn -1:n suuruinen. Pienemmillה
% luvuilla tapahtuu jokin pyצristysvirhe.

if yks>=0
    digit = rem(num, 10^(yks+1));
    digit = floor(digit/(10^yks));
else
    digit = num*10;
    digit = floor(rem(digit,10));
end
digit = num2str(digit);


function mjono = kldiv2str(div)
mjono = '      ';
if abs(div)<100
    %Ei tarvita e-muotoa
    mjono(6) = num2str(rem(floor(div*1000),10));
    mjono(5) = num2str(rem(floor(div*100),10));
    mjono(4) = num2str(rem(floor(div*10),10));
    mjono(3) = '.';
    mjono(2) = num2str(rem(floor(div),10));
    arvo = rem(floor(div/10),10);
    if arvo>0
        mjono(1) = num2str(arvo);
    end

else
    suurinYks = floor(log10(div));
    mjono(6) = num2str(suurinYks);
    mjono(5) = 'e';
    mjono(4) = palautaYks(abs(div),suurinYks-1);
    mjono(3) = '.';
    mjono(2) = palautaYks(abs(div),suurinYks);
end


%--------------------------------------------------------------------------
%Seuraavat kolme funktiota liittyvat alkupartition muodostamiseen.

function initial_partition=admixture_initialization(data_matrix,nclusters, Z)

size_data=size(data_matrix);
nloci=size_data(2)-1;
n=max(data_matrix(:,end));
T=cluster_own(Z,nclusters);
initial_partition=zeros(size_data(1),1);
for i=1:n
    kori=T(i);
    here=find(data_matrix(:,end)==i);
    for j=1:length(here)
        initial_partition(here(j),1)=kori;
    end
end

function T = cluster_own(Z,nclust)
true=logical(1);
false=logical(0);
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

function [sumcounts, counts, logml] = ...
    initialCounts(partition, data, npops, rows, noalle, adjprior)

nloci=size(data,2);
ninds = size(rows, 1);

koot = rows(:,1) - rows(:,2) + 1;
maxSize = max(koot);

counts = zeros(max(noalle),nloci,npops);
sumcounts = zeros(npops,nloci);
for i=1:npops
    for j=1:nloci
        havainnotLokuksessa = find(partition==i & data(:,j)>=0);
        sumcounts(i,j) = length(havainnotLokuksessa);
        for k=1:noalle(j)
            alleleCode = k;
            N_ijk = length(find(data(havainnotLokuksessa,j)==alleleCode));
            counts(k,j,i) = N_ijk;
        end
    end
end

%initializeGammaln(ninds, maxSize, max(noalle));

logml = laskeLoggis(counts, sumcounts, adjprior);

%--------------------------------------------------------------------------


function [partitionSummary, added] = addToSummary(logml, partitionSummary, worstIndex)
% Tiedetההn, ettה annettu logml on isompi kuin huonoin arvo
% partitionSummary taulukossa. Jos partitionSummary:ssה ei vielה ole
% annettua logml arvoa, niin lisהtההn worstIndex:in kohtaan uusi logml ja
% nykyistה partitiota vastaava nclusters:in arvo. Muutoin ei tehdה mitההn.

apu = find(abs(partitionSummary(:,2)-logml)<1e-5);
if isempty(apu)
    % Nyt lצydetty partitio ei ole vielה kirjattuna summaryyn.
    global PARTITION;
    npops = length(unique(PARTITION));
    partitionSummary(worstIndex,1) = npops;
    partitionSummary(worstIndex,2) = logml;
    added = 1;
else
    added = 0;
end

%--------------------------------------------------------------------------

function inds = returnInOrder(inds, pop, globalRows, data, ...
    adjprior, priorTerm)
% Palauttaa yksilצt jהrjestyksessה siten, ettה ensimmהisenה on
% se, jonka poistaminen populaatiosta pop nostaisi logml:n
% arvoa eniten.

global COUNTS;      global SUMCOUNTS;
ninds = length(inds);
apuTaulu = [inds, zeros(ninds,1)];

for i=1:ninds
    ind =inds(i);
    rows = globalRows(i,1):globalRows(i,2);
    diffInCounts = computeDiffInCounts(rows, size(COUNTS,1), size(COUNTS,2), data);
    diffInSumCounts = sum(diffInCounts);

    COUNTS(:,:,pop) = COUNTS(:,:,pop)-diffInCounts;
    SUMCOUNTS(pop,:) = SUMCOUNTS(pop,:)-diffInSumCounts;
    apuTaulu(i, 2) = computePopulationLogml(pop, adjprior, priorTerm);
    COUNTS(:,:,pop) = COUNTS(:,:,pop)+diffInCounts;
    SUMCOUNTS(pop,:) = SUMCOUNTS(pop,:)+diffInSumCounts;
end
apuTaulu = sortrows(apuTaulu,2);
inds = apuTaulu(ninds:-1:1,1);

%--------------------------------------------------------------------------

function [emptyPop, pops] = findEmptyPop(npops)
% Palauttaa ensimmהisen tyhjהn populaation indeksin. Jos tyhjiה
% populaatioita ei ole, palauttaa -1:n.

global PARTITION;
pops = unique(PARTITION)';
if (length(pops) ==npops)
    emptyPop = -1;
else
    popDiff = diff([0 pops npops+1]);
    emptyPop = min(find(popDiff > 1));
end
