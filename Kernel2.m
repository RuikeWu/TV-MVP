% This function is the kernel function with the form 3/4 (1-x^2)*1(abs(x) < 1)
% using Hong 2005 to eliminate boundary effect
function result = Kernel2(T,t,x,h)
    d1 = (t - x)/(T*h); 
    f= @(d) 0.75.*(1-d.*d).*(abs(d)<=1);
	if x < T*h
		int1 = integral(f,-x/(T*h),1);
		result = f(d1)/int1;
	else
		if x > (T- T*h)
		int1 = integral(f, -1, (1-x/T)/h);
		result = f(d1)/int1;
		else
		result = f(d1);
		end
	end
    
end