function covint = covintel_0(n,p,eigensamples,r1,theta1,r2,theta2)
c=p/n;
z1=r1.*exp(1i.*theta1);
z2=r2.*exp(1i.*theta2);
%%%%%%%%%%%%%%%%%calculations of sz1,sz2 and derisz1,derisz2%%%%%%%%%%%%%%%
stieltjes1=0;
stieltjes2=0;
deristieltjes1=0;
deristieltjes2=0;
for j=1:p   %calculate m(z) 
    stieltjes1=stieltjes1+1./(eigensamples(j)-z1);
    stieltjes2=stieltjes2+1./(eigensamples(j)-z2);
    deristieltjes1=deristieltjes1+1./((eigensamples(j)-z1).^2);
    deristieltjes2=deristieltjes2+1./((eigensamples(j)-z2).^2);
end
stieltjes1=stieltjes1./p;
stieltjes2=stieltjes2./p;
ulstj1=-(1-c)./z1+c.*stieltjes1; %calculate mbar(z)
ulstj2=-(1-c)./z2+c.*stieltjes2;
deristieltjes1=deristieltjes1./p;  % m(z)^' ' means derivative
deristieltjes2=deristieltjes2./p;
deulstj1=(1-c)./z1.^2+c.*deristieltjes1; %mbar(z)^'
deulstj2=(1-c)./z2.^2+c.*deristieltjes2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%limitvar%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
covint=-z1.*z2.*(4.*deulstj1.*deulstj2./(ulstj1-ulstj2).^2).*z1.*z2;
end
