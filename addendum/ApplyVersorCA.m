function tp = ApplyVersorCA(T,ca)
% TP = APPLYVERSORCA(T,CA)
%  Apply a versor to all the elements of a cell array

	for i=1:size(ca,2)
		tp{i} = zeroepsilons(T*ca{i}*inverse(T));
	end
