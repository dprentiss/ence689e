subplot(3,2,1)
plot(RFE,conceptions)
axis('tight')
title('a) Conceptions/AF')

subplot(3,2,2)
plot(RFE,salesMale)
axis('tight')
title('b) Sales/Male')

subplot(3,2,3)
plot(RFE,salesFemale)
axis('tight')
title('c) Sales/Female')

subplot(3,2,4)
plot(RFE,mortMat)
axis('tight')
title('d) Deaths/Mature Idv.')

subplot(3,2,5)
plot(RFE,mortImm)
axis('tight')
title('e) Deaths/Immature Idv.')
