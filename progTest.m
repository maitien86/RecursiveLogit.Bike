global Op;
Op.x = [-2.6,-1,-0.3,-4]';
[LL,grad,H,Hs] = getAnalyticalHessian_test();
H
[val1, gr1] = getLL_test();
h = 0.000000001 * [1 0 0 0]';
Op.x = Op.x +h;
[val2, gr2] = getLL_test();
(gr2 - gr1)/0.000000001
