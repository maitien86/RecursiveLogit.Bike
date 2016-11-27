
% get 1/mu * Grad(V0)
%%
function [gradV0] = getGradV0(M, Att, op, expV, origin)
    global isLinkSizeInclusive;
    gradExpV = getGradExpV(M, Att, op, expV);
    gradV0 = zeros(1,op.n);
    if (isLinkSizeInclusive == true)
        for i = 1:op.n
            gradV0(i) =  gradExpV(i).value(origin)/expV(origin); %write .value for link size and .Value else
        end
    else
        for i = 1:op.n
            gradV0(i) =  gradExpV(i).Value(origin)/expV(origin); %write .value for link size and .Value else
        end
    end
end

function [gradExpV] = getGradExpV(M, Att, op, expV)
    I = speye(size(M));  
    A = I - M; 
    global isLinkSizeInclusive;
    if (isLinkSizeInclusive == true)       
        gradExpV = objArray(op.n); %I added that for link size(maelle)
        for i = 1:op.n
            u = M .* (Att(i).value); %write .value for link size and .Value else
            v = u * expV; 
            gradExpV(i).value = A\v; %write gradExpV(i).value = A\v for link size and gradExpV(i) = Matrix2D(A\v) else;
        end
    else
        for i = 1:op.n
            u = M .* (Att(i).Value); %write .value for link size and .Value else
            v = u * expV; 
            gradExpV(i) = Matrix2D(A\v); %write gradExpV(i).value = A\v for link size and gradExpV(i) = Matrix2D(A\v) else;
        end
    end
end













































































