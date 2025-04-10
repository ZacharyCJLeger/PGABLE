function dtri = cgadelaunay(pts,tri)
ntris = size(tri,1);
nvts = size(pts,2)

GAScene.usefigure();
GAScene.view(2);


function pltpts(pts)
  pnpts = size(pts,2);
  for t=1:pnpts
    x(t) = pts(t).getx();
    y(t) = pts(t).gety();
    z(t) = pts(t).getz();
    text(x(t)+0.02,y(t),num2str(t));
  end
  plot3(x,y,z,'o');
end

function plttris(tris)
  ntris = size(tris,1);
  for t=1:ntris
    x(1) = pts(tris(t,1)).getx(); y(1) = pts(tris(t,1)).gety(); z(1) = pts(tris(t,1)).getz();
    x(2) = pts(tris(t,2)).getx(); y(2) = pts(tris(t,2)).gety(); z(2) = pts(tris(t,2)).getz();
    x(3) = pts(tris(t,3)).getx(); y(3) = pts(tris(t,3)).gety(); z(3) = pts(tris(t,3)).getz();
    x(4) = x(1); y(4) = y(1); z(4) = z(1);
    plot3(x,y,z);
  end
end

function pltBldTris(tris)
  ntris = size(tris,1);
  for t=1:ntris
    x(1) = pts(tris(t,1)).getx(); y(1) = pts(tris(t,1)).gety(); z(1) = pts(tris(t,1)).getz();
    x(2) = pts(tris(t,2)).getx(); y(2) = pts(tris(t,2)).gety(); z(2) = pts(tris(t,2)).getz();
    x(3) = pts(tris(t,3)).getx(); y(3) = pts(tris(t,3)).gety(); z(3) = pts(tris(t,3)).getz();
    x(4) = x(1); y(4) = y(1); z(4) = z(1);
    plot3(x,y,z,'LineWidth',3);
  end
end

function em = computeEdgeTriMatrix(tris,nv)
  nt = size(tris,1);
  em = zeros(nv,nv);
  for t=1:nt
    em(tris(t,1),tris(t,2)) = t;
    em(tris(t,2),tris(t,3)) = t;
    em(tris(t,3),tris(t,1)) = t;
  end
end

  pclf; 
  hold on;
  view(2); GAScene.view(0); axis([-2 3 -2 3]); axis equal;
  pltpts(pts);
  plttris(tri);  
  
  pause(2);
done = 0;
while ~done
  done=1;
  etm = computeEdgeTriMatrix(tri,nvts);
  for i=1:nvts
    for j=1:i-1
      if etm(i,j)==0 || etm(j,i)==0
        continue;
      end
      fprintf("Edge %d,%d\n",i,j);
      tri1 = etm(i,j); tri2 = etm(j,i);
      pt11 = tri(tri1,1); pt12 = tri(tri1,2); pt13 = tri(tri1,3);
      pt21 = tri(tri2,1); pt22 = tri(tri2,2); pt23 = tri(tri2,3);
      fprintf("%d %d: t1 %d %d %d, t2 %d %d %d\n",i,j,pt11,pt12,pt13, pt21,pt22,pt23);
  
      if i~=pt11 && j~=pt11
        ep1 = pt11; ep1i = 1; prev1 = 3;
      elseif i~=pt12 && j~=pt12
        ep1 = pt12; ep1i = 2; prev1 = 1;
      else
        ep1 = pt13; ep1i = 3; prev1 = 2;
      end
      
      if i~=pt21 && j~=pt21
        %fprintf("Case 1\n");
        ep2 = pt21; ep2i = 1; prev2 = 3;
      elseif i~=pt22 && j~=pt22
        %fprintf("Case 2\n");
        ep2 = pt22; ep2i = 2; prev2 = 1;
      else
        %fprintf("Case 3\n");
        ep2 = pt23; ep2i = 3; prev2 = 2;
      end
      c = pts(pt11)^pts(pt12)^pts(pt13);
      e = (pts(ep2)^c).*(no^e1^e2^ni);
  
      pclf; hold on; pltpts(pts); plttris(tri); 
      pltBldTris(tri(tri1,:));  
      pltBldTris(tri(tri2,:));  
      draw(c); 
      draw(pts(pt11)); draw(pts(pt12)); draw(pts(pt13));
      draw(pts(ep2),'g')
      pause(1)
       
      if double(e) < 0
	  fprintf(" Swapping edges %d,%d and %d,%d\n",i,j,ep1,ep2);
         % Need to "swap" common edge
         %    ep1           ep1
         %   /   \         / | \
         % pi-----pj  => pi  |  pj
         %   \   /         \ | / 
         %    ep2           ep2
         % in top triangle, prev1=pj; in bottom, prev2=pi
         %  So in top triangle, replace prev1 with ep2
         %  and in bot triangle, replace prev2 with ep1
	 if i==5 && j==1
           tri
	 end
         tri(tri1,prev1) = ep2;
         tri(tri2,prev2) = ep1;
	 done = 0;
	draw(pts(ep2),'r');  
	PGABLEDraw.plotline({pts(i),pts(j)},'Color','w','LineStyle',':','LineWidth',5);
	pause(2)
	PGABLEDraw.plotline({pts(ep1),pts(ep2)},'LineWidth',4);
	pause(1)
	pclf; hold on;pltpts(pts);plttris(tri); 
	pause(1)
      end
      if done==0
        break;
      end
    end
    if done==0
      break;
    end
  end
  pclf; hold on; pltpts(pts);plttris(tri); 
end  
dtri =tri;
end
