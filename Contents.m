% PGABLE is a Matlab toolkit for geometric algebra.
%  
%=====================================================================
%
%GA is a parent class of all geometric algebra models.
%OGA and PGA are currently implemented.
%To switch between the two models, run GA.model(OGA) and GA.model(PGA)
%Run "GA.settings", "PGA.settings" and "OGA.settings" for more settings.
%
%The method draw allows GA elements to be drawn to a figure. To see how to manage this
%figure, run "help GAScene".
% 
%Below is a summary of the operations PGA provides:
%
%PGA  is a child class of GA for elements of Projective/Plane-based Geometric Algebra.
%   Elements
%      Basic elements include e0, e1, e2, e3, e01, e02, e03, e12, e31, e23, e021,
%      e013, e032, e123, e0123.
%      Additionally, we have e13 = -e31, e012 = -e21, e023 = -e032.
%      We also have method for creating PGA points, point(x, y, z), which creates a
%      PGA point with coordinates (x, y, z). We also have origin() = point(0, 0, 0).
%
%   Operations
%      You can use these special characters for these basic operations:
%         • +  for addition               also: plus(A, B)
%         • -  for subtraction            also: minus(A, B)
%         • *  for the geometric product  also: product(A, B)
%         • /  for division               also: divide(A, B)
%         • ^  for the outer product      also: outer(A, B)
%         • .* for the inner product      also: inner(A, B)
%         • == for equality               also: eq(A, B)
%         • ~= for inequality             also: neq(A, B)
%      Additonally, there are basic operations:
%         • meet(A, B)                    to compute the meet of two multivectors
%         • join(A, B)                    to compute the join of two multivectors
%         • dual(A)                       to compute the dual
%         • inverse(A)                    to compute the inverse
%         • gradeinvolution(A)            to compute the grade involution
%         • conjugate(A)                  to compute the conjugate
%         • reverse(A)                    to compute the reverse
%         • norm(A)                       to compute the norm
%         • normalize(A)                  to normalize the multivector
%         • poincaredual(A)               to compute the poincare dual
%         • hodgedual(A)                  to compute the hodge dual
%         • inversehodgedual(A)           to compute the inverse hodge dual
%         • getx(A)                       to get the x coordinate of a PGA point
%         • gety(A)                       to get the y coordinate of a PGA point
%         • getz(A)                       to get the z coordinate of a PGA point
%         • zeroepsilons(A)               to zero-out epsilons (small errors)
%         • draw(A)                       to draw the multivector
%         (See also GAScene for more information on draw calls)
%         • grade(A, g)                   to select the grade-g component of a multivector
%         • isgrade(A, g)                 to determine if a multivector is of grade g
%      There are also more advanced operations:
%         • sqrt(A)                       to compute the square root
%         • glog(A)                       to compute the geometric log
%         • gexp(A)                       to compute the geometric exponent