function y = numerov(y0, y0Prime, h, V)
%求解初值问题的Numerov算法的程序
%运动方程为：y''=V(x)y
%要求U(x)和V(x)是关于x=0点对称的。
%V(x)只取x=>的部分
%y0     ->y(x=0)
%y0Prime->y'(x=0)
%h是步长
%V=V(x=>0)
nu= numel(V);
y = zeros(1,nu);
f = zeros(1,nu);
U = zeros(1,nu);
y(1) = y0;
% (*f(n)=u(n)+V(n)y(n)*)
f(1) = U(1) + V(1)*y(1);
a11 = 1.0 - V(2)*h*h/4.0;
a12 = V(3)*h*h/24.0;
a21 = -2.0 - 5.0*V(2)*h*h/6.0;
a22 = 1 - V(3)*h*h/12.0;
b1 = y0 + h*y0Prime + h*h*(7.0*f(1) + 6.0*U(2) - U(3))/24.0;
b2 = -y0 + h*h*(f(1) + 10.0*U(2) + U(3))/12.0;
y(2) = (a22*b1 - a12*b2)/(a11*a22 - a12*a21);
for i = 2:nu-1
    f(i) = U(i) + V(i)*y(i);
    temp1 =2.0*y(i)-y(i-1)+h*h*(U(i + 1)+10.0*f(i)+f(i-1))/12.0;
    temp2 =1.0 - V(i + 1)*h*h/12.0;
    y(i+1)=temp1/temp2;
end