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
