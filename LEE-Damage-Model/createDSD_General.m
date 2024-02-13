function DSD_PDF = createDSD_General(phi_d, I, N, q, A, p)

n = N.*I.^q;
a = A.*I.^p;

DSD_PDF = (n./a).*(phi_d./a).^(n-1).*exp(-(phi_d./a).^n);

end