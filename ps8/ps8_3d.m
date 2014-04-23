digits(4);
Cyy = cov(squeeze(Y_open(:,:,97,:))')
latex(vpa(sym(Cyy)))
Corryy = corr(squeeze(Y_open(:,:,97,:))')
latex(vpa(sym(Corryy)))

Cy0_pix = [1^2 0 0 0; 0 1^2 0 0; 0 0 0.05^2 0; 0 0 0 0.05^2];
for i = 1:4
    for j = 1:4
        yic(i,j) = Corryy(i,j)*sqrt(Cy0_pix(i,i))*sqrt(Cy0_pix(j,j));
    end
end

Cy0_pix = yic;

save ic3d Cy0_pix;