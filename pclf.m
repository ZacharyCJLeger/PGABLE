function pclf
	%PCLF - Clear all the scene items including vanishing objects.
	[a b] = view();
	GAScene.clearitems()
	clf;
	view([a b]);

