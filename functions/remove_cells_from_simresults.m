% remove_cells_from_simresults.m
function newsimresults = remove_cells_from_simresults(simresults, selected)

newsimresults = simresults;

newsimresults.networksize = [1 length(selected) 1];

newsimresults.networkParameters.connectivityMatrix  = ones(length(selected))*eps;

for f = fields(simresults.networkHistory)'
	str = ['newsimresults.networkHistory.' f{1} ' = newsimresults.networkHistory.' f{1} '(selected,:)';];
	eval(str)
end

for f = fields(simresults.cellParameters)'
	str = ['newsimresults.cellParameters.' f{1} ' = newsimresults.cellParameters.' f{1} '(selected,:)']
	eval(str)
end
