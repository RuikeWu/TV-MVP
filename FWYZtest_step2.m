function [Ln,ifa]=  FWYZtest_step2(X,alpha1)
    tic;
    [p,T]  = size(X);
    c=p/T;
    r1=10;
    r2=20;
    % divide sample
    X1 = X(:,1:T/2);
    X2 = X(:,T/2+1:T);
    % Calculate mean part
    eigenvecRnmean1=sort(eig((2/T).*X1*X1'));
    eigenvecRnmean2=sort(eig((2/T).*X2*X2'));
    xpart =sum(eigenvecRnmean1)-sum(eigenvecRnmean2);
    x2part =sum(eigenvecRnmean1.^2)-sum(eigenvecRnmean2.^2);

    % Calculate COV
    S1 = (2/T)*(X1*X1');
    eigensamples1 = eig(S1);
    limitvar0=-1./(4.*(pi.^2)).*integral2(@(theta1,theta2)covintel_0(T,p,eigensamples1,r1,theta1,r2,theta2),0,2*pi,0,2*pi);
    S2 = (2/T)*(X2*X2');
    eigensamples2 = eig(S2);
    limitvar2=-1./(4.*(pi.^2)).*integral2(@(theta1,theta2)covintel_2(T,p,eigensamples2,r1,theta1,r2,theta2),0,2*pi,0,2*pi);
    limitvar1=-1./(4.*(pi.^2)).*integral2(@(theta1,theta2)covintel(T,p,eigensamples1,eigensamples2,r1,theta1,r2,theta2),0,2*pi,0,2*pi);
    tilde_omega = 2.*[limitvar0,limitvar1;limitvar1,limitvar2];
    % Compute Ln
    Ln = [xpart,x2part]* inv(real(tilde_omega))*[xpart;x2part];
    chivalue =chi2inv((1-alpha1),2);
    ifa = (Ln>chivalue);
    toc;
end


