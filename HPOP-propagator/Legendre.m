% fi [rad]
function [pnm, dpnm] = Legendre(n,m,fi)

pnm = zeros(n+1,m+1);
dpnm = zeros(n+1,m+1);

pnm(1,1)=1;
dpnm(1,1)=0;
pnm(2,2)=sqrt(3)*cos(fi);
dpnm(2,2)=-sqrt(3)*sin(fi);
% diagonal coefficients
for i=2:n    
    pnm(i+1,i+1)= sqrt((2*i+1)/(2*i))*cos(fi)*pnm(i,i);
end
for i=2:n
    dpnm(i+1,i+1)= sqrt((2*i+1)/(2*i))*((cos(fi)*dpnm(i,i))- ...
                  (sin(fi)*pnm(i,i)));
end
% horizontal first step coefficients
for i=1:n
    pnm(i+1,i)= sqrt(2*i+1)*sin(fi)*pnm(i,i);
end
for i=1:n
    dpnm(i+1,i)= sqrt(2*i+1)*((cos(fi)*pnm(i,i))+(sin(fi)*dpnm(i,i)));
end
% horizontal second step coefficients
j=0;
k=2;
while(1)
    for i=k:n        
        pnm(i+1,j+1)=sqrt((2*i+1)/((i-j)*(i+j)))*((sqrt(2*i-1)*sin(fi)*pnm(i,j+1))...
                    -(sqrt(((i+j-1)*(i-j-1))/(2*i-3))*pnm(i-1,j+1)));
    end
    j = j+1;
    k = k+1;
    if (j>m)
        break
    end
end
j = 0;
k = 2;
while(1)
    for i=k:n        
        dpnm(i+1,j+1)=sqrt((2*i+1)/((i-j)*(i+j)))*((sqrt(2*i-1)*sin(fi)*dpnm(i,j+1))...
             +(sqrt(2*i-1)*cos(fi)*pnm(i,j+1))-(sqrt(((i+j-1)*(i-j-1))/(2*i-3))*dpnm(i-1,j+1)));
    end
    j = j+1;
    k = k+1;
    if (j>m)
        break
    end
end

