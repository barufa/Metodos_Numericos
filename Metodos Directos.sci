function [L,U] = doolitle(A,n)
    L=zeros(n,n)
    U=zeros(n,n)
    for i=1:n
        for j=1:n
            if i<=j then
                U(i,j)=A(i,j)
                for k=1:i-1
                    U(i,j) = U(i,j) - L(i,k)*U(k,j);
                end
            end
            if j<=i then
                L(i,j) = A(i,j);
                for k=1:j-1
                    L(i,j) = L(i,j) - L(i,k)*U(k,j)
                end
                L(i,j) = L(i,j)/U(j,j);
            end
        end
    end
endfunction

function U = cholesky(A,n)
    U=zeros(n,n)
    for i=1:n
        for j=i:n
            if i==j then
                t=0
                for k=1:(i-1)
                    t=t+(U(k,i)*U(k,i))
                end
                U(j,j)=sqrt(A(i,i)-t)
            else 
                t=-A(i,j)
                for k=1:i-1
                    t=t+(U(k,i)*U(k,j))
                end
                U(i,j)=-t/U(i,i)
            end
        end
    end
endfunction

function S = sistema(A,n,v)
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

function S = resolver_dolittle(A,n,v)
    [L,U] = doolitle(A,n)
    v=inv(L)*v'
    S=sistema(U,n,v)
endfunction

function S = resolver_cholesky(A,n,v)
    U = cholesky(A,n)
    v=inv(U')*v'
    S=sistema(U,n,v)
endfunction

