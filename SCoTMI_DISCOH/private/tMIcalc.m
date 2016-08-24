function tMI = tMIcalc(partCoh, NFFT, Tm)
%%
tMIscaling = -Tm./(2*NFFT);
tMIbuf = tMIscaling * log(1 - abs(partCoh).^2);
tMI = sum(tMIbuf, 3);


end