function hilbertBspline(order,cpf)
% Hilbert generation from https://blogs.mathworks.com/steve/2012/01/25/generating-hilbert-curves/
if ~exist('cpf','var')
  cpf=0;
end

A = zeros(0,2);
B = zeros(0,2);
C = zeros(0,2);
D = zeros(0,2);

north = [ 0  1];
east  = [ 1  0];
south = [ 0 -1];
west  = [-1  0];

%order = 5;
for n = 1:order
  AA = [B ; north ; A ; east  ; A ; south ; C];
  BB = [A ; east  ; B ; north ; B ; west  ; D];
  CC = [D ; west  ; C ; south ; C ; east  ; A];
  DD = [C ; south ; D ; west  ; D ; north ; B];

  A = AA;
  B = BB;
  C = CC;
  D = DD;
end

A = [0 0; cumsum(A)];

[r,c] = size(A);
for i=1:r
  hcp(i) = gapoint(A(i,1),A(i,2),0);
end
drawbspline(2,hcp,(1:r+2),cpf)
