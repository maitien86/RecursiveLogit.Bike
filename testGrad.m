
global Op;
step = 0.0000001;
%Op.x = [-1;-1;-1;-2;-0.2;0.3;0.5];
%Op.x = [-2.4,-1,-0.4,-4.4, -0.5, -0.7]';
x = Op.x;
[val1, gr1] = getLL_test();
gr1
H  = eye(Op.n);
for i = 1:4
    h = step * H(:,i);
%h = step * [0 0 0 0 0 1]';
    Op.x = x +h;
    [val2, gr2] = getLL_test();
    (val2 - val1)/step
end
