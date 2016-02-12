B=reshape(MULT,1,size(MULT,1)*size(MULT,2));
BB=abs(B(~isnan(B)));

norm(BB)/sqrt(length(BB))
max(BB)
median(BB)
100*sum(BB<0.1)/length(BB)
100*sum(BB>=0.1 & BB<0.3)/length(BB)
100*sum(BB>=0.3 & BB<0.5)/length(BB)
100*sum(BB>=0.5 & BB<1)/length(BB)
100*sum(BB>=1)/length(BB)