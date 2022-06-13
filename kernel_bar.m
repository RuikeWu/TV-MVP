% This function is the kernel bar function of Epanechnikov kernel function
function result = kernel_bar(T,t,x,h)
    d = (t - x)/(T*h);
    result = (3/5 - (3/4) * d^2 + (3/8)*(abs(d)^3)-(3/160)*(abs(d)^5))*(abs(d)<=2);
end