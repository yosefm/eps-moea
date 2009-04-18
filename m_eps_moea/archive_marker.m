% Selects for breeding an individual from the archive.
% Currently selects randomly.
%
% Arguments:
% archive_marker - a length-p boolean vector stating which of the population
%	is in the archive.
%
% Returns:
% sel - the index in the population of the selected archive member.

function sel = archive_select(archive_marker)
	archive_size = sum(archive_marker);
	archive_ind = find(archive_marker);
	sel = archive_ind(fix(rand() * archive_size + 1));
end
