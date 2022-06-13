% R is N by T
% F is m by T 
function  [obj_sum] = V(R,F,Lambda,FF)
 %       obj_sum = 0;
        [N,T] = size(R);
        m = size(F,1);
        %FF = zeros((m*T),T);
        for i = 1 :T 
            FF((m*i - (m-1)): (m*i),i) = F(:,i);
        end
        BB = [];
        for r = 1 : T
           BB = [BB,(Lambda{r})];
        end
        obj_sum =  norm(R - BB*FF,'fro');
        
%         for i = 1 : size(R,1)
%             Lambda_set = [];
%             for r = 1 :size(R,2);
%                Lambda_temp = Lambda{r};
%                Lambda_set = [Lambda_set,Lambda_temp(i,:)']; 
%             end
%             obj_sum = obj_sum +sum(diag(diag(R(i,:))- F* Lambda_set).^2);
% %             for r = 1: size(R,2)
% %                 Lambda_temp = Lambda{r};
% %                 obj_sum = obj_sum + (R(i,r) - F(r,:)*Lambda_temp(i,:)')^2;
% %             end
%         end
end
    