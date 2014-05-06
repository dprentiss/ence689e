%subplot(3,2,1)
figure
plot(RFE,conceptions)
axis('tight')
title('Conceptions/AF')
xlabel('RFE')

%subplot(3,2,2)
figure
plot(RFE,salesMale)
axis('tight')
title('Sales/Male')
xlabel('RFE')

%subplot(3,2,3)
figure
plot(RFE,salesFemale)
axis('tight')
title('Sales/Female')
xlabel('RFE')

%subplot(3,2,4)
figure
plot(RFE,mortMat)
axis('tight')
title('Deaths/Mature Idv.')
xlabel('RFE')

%subplot(3,2,5)
figure
plot(RFE,mortImm)
axis('tight')
title('Deaths/Immature Idv.')
xlabel('RFE')
