function x = sumar_fila_tri(a,i,mult,n)
    for k=i:i+2
        if k==n+1 then
            break
        end
        a(i+1,k)=a(i+1,k)-(mult*a(i,k))
    end
    x=a
endfunction

function [l,a] = gauss_tridiagonal(a,n)
    l = eye(n,n)
    for i=1:n-1
        l(i+1,i) = a(i+1,i)/a(i,i)
        a=sumar_fila_tri(a,i,l(i+1,i),n)
    end
endfunction

function S = slove_tridiagonal(A,n,vcons)
    S = zeros(n,1)
    if det(A) then
        [L,A] = gauss_tridiagonal(A,n)
        vcons = inv(L) * vcons'
        S(n)=vcons(n)/A(n,n)
        for i = 1:n-1
            s=0
            e=n-i
            for j= e+1:n
    
            s=s+A(e,j)*S(j)
            end
            S(e)=(vcons(e)-s)/A(e,e)
        end
    else
        printf("Error: --Matriz singular\n")
    end
endfunction
