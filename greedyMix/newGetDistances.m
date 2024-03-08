function [Z, dist] = newGetDistances(data, rowsFromInd)

  ninds = max(data(:,end));
  nloci = size(data,2)-1;
  riviLkm = nchoosek(ninds,2);

  empties = find(data<0);
  data(empties)=0;
  data = uint8(data);   % max(noalle) oltava <256

  pariTaulu = zeros(riviLkm,2);
  aPointer=1;
  for a=1:ninds-1
    pariTaulu(aPointer:aPointer+ninds-1-a,1) = ones(ninds-a,1)*a;
    pariTaulu(aPointer:aPointer+ninds-1-a,2) = (a+1:ninds)';
    aPointer = aPointer+ninds-a;
  end

  eka = pariTaulu(:,ones(1,rowsFromInd));
  eka = eka * rowsFromInd;
  miinus = repmat(rowsFromInd-1 : -1 : 0, [riviLkm 1]);
  eka = eka - miinus;

  toka = pariTaulu(:,ones(1,rowsFromInd)*2);
  toka = toka * rowsFromInd;
  toka = toka - miinus;

  %eka = uint16(eka);
  %toka = uint16(toka);

  summa = zeros(riviLkm,1);
  vertailuja = zeros(riviLkm,1);

  clear pariTaulu; clear miinus;

  x = zeros(size(eka));    x = uint8(x);
  y = zeros(size(toka));   y = uint8(y);

  for j=1:nloci;

    for k=1:rowsFromInd
      x(:,k) = data(eka(:,k),j);
      y(:,k) = data(toka(:,k),j);
    end

    for a=1:rowsFromInd
      for b=1:rowsFromInd
        vertailutNyt = double(x(:,a)>0 & y(:,b)>0);
        vertailuja = vertailuja + vertailutNyt;
        lisays = (x(:,a)~=y(:,b) & vertailutNyt);
        summa = summa+double(lisays);
      end
    end
  end

  clear x;    clear y;   clear vertailutNyt;
  nollat = find(vertailuja==0);
  dist = zeros(length(vertailuja),1);
  dist(nollat) = 1;
  muut = find(vertailuja>0);
  dist(muut) = summa(muut)./vertailuja(muut);
  clear summa; clear vertailuja;

  Z = linkage(dist');


end
