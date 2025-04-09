function tri = cgatriangulate(pts)
npts = size(pts,2);

function pltpts(pts)
  pnpts = size(pts,2);
  for t=1:pnpts
    x(t) = pts(t).getx();
    y(t) = pts(t).gety();
    z(t) = pts(t).getz();
    text(x(t)+0.03,y(t),num2str(t));
  end
  plot3(x,y,z,'o');
end

function plttris(tris)
  ntris = size(tris,1);
  for t=1:ntris
    %fprintf("tri %d %d %d\n",tris(t,1),tris(t,2),tris(t,3));
    x(1) = pts(tris(t,1)).getx(); y(1) = pts(tris(t,1)).gety(); z(1) = pts(tris(t,1)).getz();
    x(2) = pts(tris(t,2)).getx(); y(2) = pts(tris(t,2)).gety(); z(2) = pts(tris(t,2)).getz();
    x(3) = pts(tris(t,3)).getx(); y(3) = pts(tris(t,3)).gety(); z(3) = pts(tris(t,3)).getz();
    x(4) = x(1); y(4) = y(1); z(4) = z(1);
    plot3(x,y,z);
  end
end

function em = computeEdgeMatrix(tris,nv)
  nt = size(tris,1);
  em = zeros(nv,nv);
  for t=1:nt
    em(tris(t,1),tris(t,2)) = 1;
    em(tris(t,2),tris(t,3)) = 1;
    em(tris(t,3),tris(t,1)) = 1;
  end
end

function [v1a,v2b] = getNextPrevEdge(em,v1,v2)
[r,c] = size(em);
    for j=1:r
      if j~=v2 && (em(j,v1)+em(v1,j)==1)
        v1a = j;
      end
      if j~=v1 && (em(j,v2)+em(v2,j)==1)
        v2b = j;
      end
    end
end

function b = pointInTriangle(P, P1,P2,P3)
  % Construct plane containing Pi,Pj,ni, and e3
  pln12 = P1^P2^ni^e3;
  pln23 = P2^P3^ni^e3;
  pln31 = P3^P1^ni^e3;
  if double((P^pln12).*I5) > 0 && ...
     double((P^pln23).*I5) > 0 && ...
     double((P^pln31).*I5) > 0
      b=1;
   else
      b=0;
   end
end

GAScene.view(2);
pclf; hold on;
pltpts(pts); view([0 90]); axis equal; hold on;
axis([-2 3 -2 3]);

tri(1,1) = 1;
tri(1,2) = 2;
tri(1,3) = 3;
ntris = 1;

c = pts(1)^pts(2)^pts(3);
%double(c.*(ni^e1^e2)) 
if double(c.*(ni^e1^e2)) < 0
  tri(1,2) = 3;
  tri(1,3) = 2;
  c = -1*c;
end

pts(tri(1,1));
pts(tri(1,2));
pts(tri(1,3));

em = computeEdgeMatrix(tri,3);

plttris(tri)

for nxtV=4:npts
  em = computeEdgeMatrix(tri,nxtV-1);
  didsplit=0; % set to true if we insert point nxtV on interior of triangle
  for j=1:ntris
    % Construct circle for triangle j
    % Check if point nxtV is inside the triangle
    if pointInTriangle(pts(nxtV), pts(tri(j,1)),pts(tri(j,2)),pts(tri(j,3)))
      % Split triangle 1:3
      fprintf("tri %d %d %d, %d\n",tri(j,1),tri(j,2),tri(j,3),nxtV);
      fprintf("Point %d in triangle %d,%d,%d; split 1:3.\n",nxtV,tri(j,1),tri(j,2),tri(j,3));
      ptA = tri(j,1); ptB = tri(j,2); ptC = tri(j,3);
      tri(j,1) = ptA; tri(j,2) = ptB; tri(j,3) = nxtV;
      ntris = ntris+1;
      tri(ntris,1) = ptB; tri(ntris,2) = ptC; tri(ntris,3) = nxtV;
      ntris = ntris+1;
      tri(ntris,1) = ptC; tri(ntris,2) = ptA; tri(ntris,3) = nxtV;
      didsplit=1;
      break;
    end
  end

  % The new vertex wasn't inside an old one; so find where to put it outside
  if ~didsplit
    fprintf("Point outside all triangles; add to exterior.\n");
    em = computeEdgeMatrix(tri,nxtV-1);
    % Iterate through edges, looking for exterior edges
    done = 0; % used when we find the edge we're looking for
    for j=1:nxtV-1
      for k=1:nxtV-1
        % Look for an exterior edge
        if em(j,k)+em(k,j)==1
	  if em(j,k)
	    v1=j; v2=k;
	  else
	    v1=k; v2=j;
	  end
	  % Test if new point on "outside" of exterior edge
	  if double((pts(nxtV)^pts(v2)^pts(v1)^ni^e3).*I5) > 0
	     ntris = ntris+1;
	     tri(ntris,1) = v2; tri(ntris,2) = v1; tri(ntris,3) = nxtV;
             done = 1;
	  end
	end
	if done
	  break;
	end
      end
      if done
         break;
      end
    end
    plttris(tri)
    % We want the exterior of the triangulation to be convex
    % So check triangles formed with next two edges to make
    %  sure everything is convex, adding a new triangle if
    %  something is concave.  Make sure we don't recompute
    %  em before this first call
    [v1a,v2b] = getNextPrevEdge(em,v1,v2);;
    c1 = pts(nxtV)^pts(v1)^pts(v1a);
    c2 = pts(v2b)^pts(v2)^pts(nxtV);
    if double(c1.*(ni^e1^e2)) > 0
      ntris = ntris+1;
      tri(ntris,1) = v1; tri(ntris,2) = v1a; tri(ntris,3) = nxtV;
    end
    if double(c2.*(ni^e1^e2)) > 0
      ntris = ntris+1;
      tri(ntris,1) = v2; tri(ntris,2) = v2b; tri(ntris,3) = nxtV;
    end
  end
plttris(tri);
pause(1)
end


end
