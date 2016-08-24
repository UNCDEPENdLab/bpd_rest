

function discohproj_raw = DISCOHprojected_calc(partCoh, NFFT)
%#codegen
% coder.inline('never')

scaling = 1./(NFFT);

discohproj_buf = scaling .* sqrt((1 - (abs(partCoh)).^2 ));

discohproj_raw = sum(discohproj_buf,3);


% tMIscaling = -Tm./(2*NFFT);
% tMIbuf = tMIscaling * log(1 - abs(partCoh).^2);
% tMI = sum(tMIbuf, 3);


end