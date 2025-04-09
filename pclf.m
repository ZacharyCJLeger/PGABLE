function pclf
	%PCLF - Clear all the scene items including vanishing objects.
 	[a b] = view();
	ax = axis;
	GAScene.clearitems()
	clf;
	axis(ax);
	view([a b]);

