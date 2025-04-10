 pts(1) = gapoint(0,0,0);
 pts(2) = gapoint(0,1,0);
 pts(3) = gapoint(1,0,0);
 pts(4) = gapoint(0.25,0.25,0);
 pts(5) = gapoint(-0.5,2,0);
 pts(6) = gapoint(0.8,1,0);
 pts(7) = gapoint(0.2,0.7,0);
 pts(8) = gapoint(1.4,1.4,0);
 pts(9) = gapoint(-0.6,-0.5,0);
 pts(9) = gapoint(-0.6,-0.5,0);
 pts(10) = gapoint(1,2.5,0);
 pts(11) = gapoint(2.5,0.,0);

pts2(1) = gapoint(0,0,0);
pts2(2) = gapoint(2,0,0);
pts2(3) = gapoint(1,0.5,0);
pts2(4) = gapoint(1,-0.5,0);

fprintf("Triangulate Points\n");
tri = cgatriangulate(pts(1:11));
fprintf("Triangulation finished.\n");
pause(4);
fprintf("Start Delaunay triangulation.\n");
cgadelaunay(pts(1:11),tri)
fprintf("Delaunay triangulation complete.\n");
