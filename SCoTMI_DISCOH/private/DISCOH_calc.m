function notscotmi_raw = DISCOH_calc(partCoh, NFFT)
%#codegen
% coder.inline('never')

scaling = 1./(NFFT*sqrt(2));

notscotmi_buf = scaling .* sqrt(2*(1 - abs(real(partCoh))));

notscotmi_raw = sum(notscotmi_buf,3);

end