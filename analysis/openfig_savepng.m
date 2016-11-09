% openfig_savepng.m


h = get(0,'children');

for i=1:length(h)
   export_fig(num2str(i), '-png')
   saveas(h(i), ['figure' num2str(i)], 'fig');
end

