function genInnerTable
    
IPA = zeros(32);

GA.model(CGA)
    be = [1,no,e1,e2,e3,ni, no^e1,no^e2,no^e3,no^ni, e1^e2,e1^e3,e1^ni, e2^e3,e2^ni, e3^ni,...
         no^e1^e2,no^e1^e3,no^e1^ni, no^e2^e3,no^e2^ni, no^e3^ni, e1^e2^e3,e1^e2^ni, e1^e3^ni, e2^e3^ni,...
	  no^e1^e2^e3, no^e1^e2^ni, no^e1^e3^ni, no^e2^e3^ni, e1^e2^e3^ni, no^e1^e2^e3^ni];

    bename = ...
    ["scal", ...
            "EO", "E1", "E2", "E3", "EI", ...
            "EO1", "EO2", "EO3", "EOI", "E12", "E13", "E1I", "E23", "E2I", "E3I", ...
            "EO12", "EO13", "EO1I", "EO23", "EO2I", "EO3I", "E123", "E12I", "E13I", "E23I", ...
            "EO123", "EO12I", "EO13I", "EO23I", "E123I", ...
            "EO123I"]

for i=1:32
    ben(i) = -1*be(i);
end

sz = 32;

function AddTerm(sgn,alpha,beta,gamma)
  if gamma==0
   return;
  end
  IPA(gamma,beta) = sgn*alpha;   
end

function gamma = reverseLookup(a)
    gamma=0;
    for ii=1:32
        if eeq(a,be(ii))
	       gamma = ii;
	       return;
       end
       if eeq(a,ben(ii))
	       gamma = -ii;
       end
    end
end
    
function ip = computeInner(a,b)
    ip = 0;
    for r=0:5
      for s=0:5
         rslt = grade(a,r)*grade(b,s);
% dot product
%         ip = ip+grade(rslt,abs(r-s));
% left contraction
%	 if s-r>=0 
%            ip = ip+grade(rslt,s-r);
%         end
% right contraction
	 if r-s>=0 
            ip = ip+grade(rslt,r-s);
         end
      end
    end
end

numz = 0;
numNonz=0;

 for i=1:32
    for j=1:32
        ip = computeInner(be(i),be(j));
	gamma = reverseLookup(ip);
	if gamma == 0
	  numz = numz+1;
	else
	  numNonz = numNonz+1;
	end
        fprintf("%d,%d: %-14s . %-14s = %-14s gamma %d beta %d\n",i,j,char(be(i)), char(be(j)), char(ip),gamma,j);
	AddTerm(sign(gamma),i,j,abs(gamma));
    end
 end
IPA

numz
numNonz

fid = fopen("ipm.txt","w");
for i=1:32
  for j=1:32
    if IPA(i,j) < 0
      ss = "-";
    else  
      ss = " ";
    end
    if  IPA(i,j) == 0
      fprintf(fid,"%-7s", "0");
    else
      fprintf(fid,"%1s%-6s ", ss,bename(abs(IPA(i,j))));
    end
  end
  fprintf(fid,";\n");
end
fclose(fid);
 end
