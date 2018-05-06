function x = descenso(A,b,x0,iter,e)
    i=1
    x=x0'
    x0=x
    v=b'-A*x
    t=(v'*v)/(v'*(A*v))
    x=x+t*v
    while i<=iter & (norm(x-x0,'inf'))>=e
        x0=x
        v=b'-A*x
        t=(v'*v)/(v'*(A*v))
        x=x+t*v
        i=i+1
    end
endfunction

function [v,r] = potencia(A,v,iter,e)
    k=1
    v=v'
    t=v
    y=A*v
    r=(y(1))/(v(1))
    v=y/norm(y,'inf')
    while k<=iter-1 & (norm(v-t,'inf'))>=e
        t=v
        y=A*v
        v=y/norm(y,'inf')
        k=k+1
    end
    r=y(1)/v(1)
endfunction

function S = resolver(A,n,v)
    S = zeros(n,1)
    S(n)=v(n)/A(n,n)
    for i = 1:n-1
        s=0
        e=n-i
        for j= e+1:n
            s=s+A(e,j)*S(j)
        end
        S(e)=(v(e)-s)/A(e,e)
    end
endfunction

function [x,r]=inverso(A,x,n,iter,e)
    [P,L,A]=gauss_escalable(A,n)
    k=1
    x=inv(L)*(P*x')
    t=x
    x=resolver(A,n,x)
    while k<=iter & (norm(x-t,'inf'))>=e
    t=x
    x=resolver(A,n,x)
    k=k+1
    end
endfunction

function [x,r] = inverso2(A,x,iter,e)
    [x,r]=potencia(inv(A),x,iter,e)
    r=(r^(-1))
endfunction
