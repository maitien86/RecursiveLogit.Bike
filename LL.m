function [f g] = LL(x)
    global Op;
    global isLinkSizeInclusive;
    Op.x = x;
    if (isLinkSizeInclusive == true)
        [f g] = getLL(); %write getLL() for link size and getLL_test() else
    else
        [f g] = getLL_test();
    end
    Op.nFev  = Op.nFev + 1;
%    Op.value = f;
%    Op.grad = g;
%    Op.value
%    Op.grad
end
