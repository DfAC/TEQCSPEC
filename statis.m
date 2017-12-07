B=reshape(SatVal_Active,1,size(SatVal_Active,1)*size(SatVal_Active,2));
BB=abs(B(~isnan(B)));

disp(['Square root mean - ' num2str(norm(BB)/sqrt(length(BB)))])
disp(['Max amplitude - ' num2str(max(BB))])
disp(['median amplitude - ' num2str(median(BB))])
disp(['<0.1m percentage - '  num2str(100*sum(BB<0.1)/length(BB))])
disp(['0.1m<x<0.3m percentage - '  num2str(100*sum(BB>=0.1 & BB<0.3)/length(BB))])
disp(['0.3m<x<0.5m percentage - '  num2str(100*sum(BB>=0.3 & BB<0.5)/length(BB))])
disp(['0.5m<x<1m percentage - '  num2str(100*sum(BB>=0.5 & BB<1)/length(BB))])
disp(['>1m percentage - '  num2str(100*sum(BB>=1)/length(BB))])

clear B BB