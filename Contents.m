%   PGABLE is a Matlab toolkit for geometric algebra.
%
%   GA is a parent class of all geometric algebra models.
%   The models OGA and PGA are currently implemented as child classes of GA.
%   There are general settings stored in GA which apply to all GA models, and there
%   are settings stored for each model of GA.
%   To access each of these settings, run "GA.settings", "OGA.settings", or
%   "PGA.settings".
%
%   ~ BASICS YOU SHOULD KNOW ABOUT ~
%   By default you are running the PGA model.
%   To switch between PGA and OGA, run GA.model(OGA) or GA.model(PGA)
%
%   To create elements of PGA or OGA, it is recommened to write equations using the basis
%   elements. For example, one could run "a = e1 + 2*e12" to create the corresponding
%   PGA element. However, if the model you are currently running is OGA, this will
%   create an OGA element of the same form. If you are unsure what model you are
%   currently using, simply run "GA.model". If you are unsure what model a particular
%   element is, run "modelname(x)" where "x" is the element you wish to investigate.
%   If you would like to consistently see the model elements belong to, consider
%   running "GA.indicate_model(true)".
%   When switching models, elements created in OGA will remain as OGA elements.
%   Similarly, elements created in PGA will remain as PGA elements. Thus, switching models
%   only changes the process of creating NEW elements. If you wish to create an OGA element
%   in PGA (perhaps because you are being lazy), you can specify which model you wish to
%   use as such "a = e1(PGA) + e12(PGA)". This will create a PGA element regardless of which
%   model you are in. Scalars automatically conform to the model of the rest of the equation.
%   Thus, we do not specify the model of scalars and "e1(PGA) + 1" will work with no issues.
%   Elements of differing models cannot be used together, i.e. if one were to have the
%   OGA element "a = e1 + e2" and the PGA element "b = e2 + e3", their sum "a + b" would
%   result in an error (after all, which model should this new element belong to?).
%   There are methods to convert OGA elements to PGA elements and vice versa, if needed.
%
%   The method "draw(x)" allows one to draw a representation of the geometric
%   interpretation of the element "x". This geometric interpretation depends on the model
%   of the element, of course, and so the result of "draw(e1 + e2)" will depend on the
%   model you are in.
%   All GA elements are drawn to a figure called "GA Scene", and this figure is managed
%   by the class GAScene. GAScene has publically accessable functions, however it
%   is generally not recommended to use them, except for:
%    - GAScene.clearitems() to clear the items of the figure
%    - GAScene.displayitems() to see a list of items in the figure
%   For more information on how to manage your GA scene, run "help GAScene".