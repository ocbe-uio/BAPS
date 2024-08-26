function pal = testaaPop(rivi)
  % pal=1, mikהli rivi alkaa jollain seuraavista
  % kirjainyhdistelmist? Pop, pop, POP. Kaikissa muissa
  % tapauksissa pal=0.

  if length(rivi)<3
    pal = 0;
    return
  end
  if (all(rivi(1:3)=='Pop') | ...
    all(rivi(1:3)=='pop') | ...
    all(rivi(1:3)=='POP'))
    pal = 1;
    return
  else
    pal = 0;
    return
  end
end
