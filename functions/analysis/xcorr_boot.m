function xc = xcorr_boot(sp1,sp2);


isi1 = diff(sp1);
isi2 = diff(sp2);

nboot = 100;

% boot method

for K = 1:nboot


	evs1 = cumsum(isi1(randperm(length(isi1))));
	evs2 = cumsum(isi2(randperm(length(isi2))));

	l = max(evs1(end), evs2(end));
	
	b1 = zeros(1,round(evs1(end))+1);
	b2 = zeros(1,round(evs2(end))+1);
	
	b1(evs1) = 1;
	b2(evs2) = 1;
	
	xc(K,:) = xcorr(b1,b2);
	
end



% gamfit method



% poiss method