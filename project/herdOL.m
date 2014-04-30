function [ herdOL ] = herdOL( numSeasons, herdInit, RFE, RFEs )
%   herdOL Propagates a herd state through time
%   ...

h = herdInit;

for i = 2:numSeasons % for every seasons
    
	% AFt+1 = AFt + YFt - sales rate(AFt) - death rate(AFt)
   	h(i,1,:) = h(i-1,1,:) + h(i-1,3,:)...
   	 - round(interp1(RFE,salesFemale,RFEs(i-1))*h(i-1,1,:))...
     - round(interp1(RFE,mortMat,RFEs(i-1))*h(i-1,1,:));

    % NBt+1 = conception rate(AFt)
    h(i,2,:) = round(interp1(RFE,conceptions,RFEs(i-1))*h(i-1,1,:));
    
    % YFt+1 = 0.5*(NBt - mortImm(NBt)) - death rate(YFt)
    h(i,3,:) = round(0.5*(h(i-1,2,:)...
     - round(interp1(RFE,mortImm,RFEs(i-1))*h(i-1,2,:))))...
     - round(interp1(RFE,mortMat,RFEs(i-1))*h(i-1,3,:));
 

    % YMt+1 = YMt + NBt - YFt+1 - sales rate(YMt) - death rate(YMt)
    h(i,4,:) = h(i-1,4,:) + h(i,3,:)...
     - round(interp1(RFE,salesMale,RFEs(i-1))*h(i-1,4,:))...
     - round(interp1(RFE,mortMat,RFEs(i-1))*h(i-1,4,:));
end

herdOL = h;

end

