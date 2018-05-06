function y = biseccion(ff,a,b,d,e,m)
    deff("y=f(x)","y="+ff)
    m=m-1
    if f(a)*f(b) < 0 then
        c = (b+a)/2
        if (b-a)<d | abs(f(c)) < e | m <= 0 | f(c) == 0 then
            y = c
        else
            if f(a)*f(c) < 0 then
                y = biseccion(ff,a,c,d,e,m-1)
            else
                y = biseccion(ff,c,b,d,e,m-1)
            end
        end
        
    else
        y = "No es posible ejecutar sobre el intervalo dado"
    end
endfunction

function [r1,r2] = newton(ff,x,d,e,m)
    deff("y=f(x)","y="+ff)
    deff("y=Df(x)","y=(numderivative(f,x))")
    a = x
    s = 0
    n = 0
    while n < m & abs(f(a)) > e & abs(a-s) > d 
        s = a
        a = a - (f(a)/Df(a))
        n = n + 1
    end
    r1 = a
    r2 = n
endfunction


function [r1,r2] = secante(ff,x0,x1,d,e,m)
    deff("y=f(x)","y="+ff)
    deff("y=Df(x)","y=(numderivative(f,x))")
    a=x1-f(x1)*(x0-x1)/(f(x0)-f(x1))
    a0=x1
    a1=a
    n=1
    while n < m & abs(f(a)) > e & abs(a0-a1) > d
        a=a1-f(a1)*(a0-a1)/(f(a0)-f(a1))
        a0=a1
        a1=a
        n=n+1
    end
    r1 = a
    r2 = n
endfunction

function y = g(x)
  f1 = x(1)^2 + x(2)^2 - 1
  f2 = x(1) - x(2)
  y = [f1; f2]
endfunction

function y = newton_no_lineal(a,d,e,m)
    n=0
    while n < m
        a = a - (numderivative(g, a))^(-1)*g(a)
        n = n + 1
    end
    y = a
endfunction

function y = puntofijo(ff,x,e,m)
    deff("y=f(x)","y="+ff)
    n = 0
    a = x
    while n < m & abs(f(a)-a) > e do
        a = f(a)
        n = n + 1
    end
    y = a
endfunction
