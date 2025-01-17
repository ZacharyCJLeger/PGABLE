ccp(1,1) = gapoint(0,0,0); % 300
ccp(2,1) = gapoint(1,0,0.5); % 210
ccp(3,1) = gapoint(2,0,0.5); % 120
ccp(4,1) = gapoint(3,0,0); % 030
ccp(1,2) = gapoint(0,1,0.5); % 201
ccp(2,2) = gapoint(1,1,0.75); % 111
ccp(3,2) = gapoint(2,1,0.5); % 021
ccp(1,3) = gapoint(0,2,0.5); % 102
ccp(2,3) = gapoint(1,2,0.5); % 012
ccp(1,4) = gapoint(0,3,0); % 003
n=3;

qcp(1,1) = (ccp(1,1)+ccp(2,1)+ccp(1,2))/3; % 200
qcp(2,1) = (ccp(2,1)+ccp(3,1)+ccp(2,2))/3; % 110
qcp(3,1) = (ccp(3,1)+ccp(4,1)+ccp(3,2))/3; % 020
qcp(1,2) = (ccp(1,2)+ccp(2,2)+ccp(1,3))/3; % 101
qcp(2,2) = (ccp(2,2)+ccp(3,2)+ccp(2,3))/3; % 011
qcp(1,3) = (ccp(1,3)+ccp(2,3)+ccp(1,4))/3; % 002
n=2;

lcp(1,1) = (qcp(1,1)+qcp(2,1)+qcp(1,2))/3; % 100
lcp(2,1) = (qcp(2,1)+qcp(3,1)+qcp(2,2))/3; % 010
lcp(1,2) = (qcp(1,2)+qcp(2,2)+qcp(1,3))/3; % 001
n=1;

n=3;
cp = ccp;

for i=1:n+1
  for j=1:n+1-i+1
    draw(cp(i,j))
  end
end

% Draw the control net
for i=1:n+1
  for j=1:n+1-i
    va = [cp(i,j), cp(i+1,j), cp(i,j+1), cp(i,j)];
    PGADrawPolyline(va,'b','LineWidth',2);
  end
end


for i=1:n
  for j=1:n+1-i
    draw(qcp(i,j),'r')
  end
end

% Draw the control net
for i=1:n
  for j=1:n+1-i-1
    va = [qcp(i,j), qcp(i+1,j), qcp(i,j+1), qcp(i,j)];
    PGADrawPolyline(va,'r','LineWidth',2);
  end
end


for i=1:n-1
  for j=1:n+1-i-1
    draw(lcp(i,j),'g')
  end
end

% Draw the control net
for i=1:n-1
  for j=1:n+1-i-2
    va = [lcp(i,j), lcp(i+1,j), lcp(i,j+1), lcp(i,j)];
    PGADrawPolyline(va,'g','LineWidth',2);
  end
end

draw((lcp(1,1)+lcp(2,1)+lcp(1,2))/3,'k');

