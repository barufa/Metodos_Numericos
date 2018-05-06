function [L,D,U] = Fact(A,n)
    D=zeros(n,n)
    L=zeros(n,n)
    U=zeros(n,n)
    for i= 1:n
        for j=1:n
            if i==j then
                D(i,j)=A(i,j)
            end
            if i<j then
                U(i,j)=A(i,j)
            end
            if i>j then
                L(i,j)=A(i,j)
            end
        end
    end
endfunction

function x = Gauss_Seidel(A,b,x0,n,iter,e)
      x=x0'
      if iter==0 then
          return
      end
      [L,D,U]=Fact(A,n)
      Q=inv(L+D)
      I=eye(n,n)
      r=norm((I-Q*A),'inf')
      if r>=1 then
          r=1
          printf("No es posible asegurar la convergencia del metodo\n");
      else
          r=r/(1-r)
      end
      i=1
      x0=x
      x=((-1*Q)*U*x)+(Q*b')
      while i<=iter & abs((norm(x-x0,'inf')*r))>=e
          x0=x
          x=((-1*Q)*U*x)+(Q*b')
          i=i+1
      end
      printf("Pasos=%d\n",i-1)
endfunction

function x = Jacobi(A,b,x0,n,iter,e)
      x=x0'
      if iter==0 then
          return
      end
      [L,D,U]=Fact(A,n)
      Q=inv(D)
      I=eye(n,n)
      r=norm((I-Q*A),'inf')
      if r>=1 then
          r=1
          printf("No es posible asegurar la convergencia del metodo\n");
      else
          r=r/(1-r)
      end
      x0=x
      i=0
      x=(I-Q*A)*x+(Q*b')
      while i<=iter & abs((norm(x-x0,'inf')*r))>=e
          x0=x
          x=(I-Q*A)*x+(Q*b')
          i=i+1
      end
      printf("Pasos=%d\n",i-1)
endfunction

function x = Richarson(A,b,x0,n,iter,e)
      x=x0'
      if iter==0 then
          return
      end
      I=eye(n,n)
      r=norm((I-A),'inf')
      if r>=1 then
          r=1
          printf("No es posible asegurar la convergencia del metodo\n");
      else
          r=r/(1-r)
      end
      x0=x
      i=1
      x=(I-A)*x+b'
      while i<=iter & abs((norm(x-x0,'inf')*r))>=e
          x0=x
          x=(I-A)*x+b'
          i=i+1
      end
      printf("Pasos=%d\n",i-1)
endfunction
