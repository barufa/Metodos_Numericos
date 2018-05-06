function print_matriz(A,n)
    for i= 1:n
        for j=1:n
        printf ("%f ",A(i,j))
        end
        printf("\n")
    end
    printf("\n")
endfunction


function x = sumar_fila(a,i,j,mult,n)
    for k=1:n
        a(i,k) = a(i,k) + (-mult)*a(j,k) 
    end
    x=a
endfunction

function y=verificar_matriz(L,n)
    y=0
    for i=1:n
        for j=1:n
            if ~(i==j) & ~(L(i,j)==0) then
                y=1
            end
        end
    end
endfunction

function a = permutacion_aux(a,n,i,j)
    for k=1:n
        t=a(i,k)
        a(i,k)=a(j,k)
        a(j,k)=t
    end
endfunction

function [a,v,l] = permutacion(a,n,j,v,l)
    k=-1
    for i = j:n
        //Existe un pivote distinto de 0 pues la matriz A no es singular
        if ~(a(i,j) == 0)//Cuando encuentra un pivote no nulo corta
            k=i
            break;
        end
    end
    //Permuta las filas de la matrices para utilizar el pivote no nulo
    a = permutacion_aux(a,n,k,j)
    v = permutacion_aux(v,n,k,j)
    l = permutacion_aux(l,n,k,j)
endfunction

function [a,p,l]= permutacion_escalable(a,n,j,p,l)
    M = -1
    t=-1
    for i=j:n
        m=-1
        for k=j:n
            if abs(a(i,k))>m then//Busca el maximo valor de la k-esima fila
                m = abs(a(i,k))
            end
        end
        if (abs(a(i,j))/m)>t then//Busca el valor que verifica la condicion del ejercicio 3
            M=i
            t=(abs(a(i,j))/m)
        end
    end
    //Permuta las filas las matrices para usar la fila M como pivote
    a=permutacion_aux(a,n,M,j)
    p=permutacion_aux(p,n,M,j)
    l=permutacion_aux(l,n,M,j)
endfunction


function [l,a] = gauss_sin_pivote(a,n)
    l = eye(n,n)
    for j=1:n-1//Aplica metodo de Gauss sin pivoteo si es posible
        if a(j,j)==0 then
            printf("Es necesario permutar filas\n")
            a=eye(n,n)
            l=eye(n,n)
        end
        for i=(j+1):n
            l(i,j) = a(i,j)/a(j,j)
            a= sumar_fila(a,i,j,l(i,j),n)
        end
    end
endfunction

function [p,l,a] = gauss_con_pivote(a,n)
    p = eye(n,n)
    l = zeros(n,n)
    if det(a) then
        for j=1:n-1
            while a(j,j) == 0//Si el pivote es nulo realiza un intercambio de filas
                [a,p,l]=permutacion(a,n,j,p,l)
            end
            for i=(j+1):n//Completa la matriz L modifica la matriz A
                l(i,j) = a(i,j)/a(j,j)
                a=sumar_fila(a,i,j,l(i,j),n)
            end
        end
        for i=1:n//Agrega los 1 en la diagonal
           l(i,i)=1
        end
    else
        printf("Matriz singular\n");
    end
endfunction

function [p,l,a] = gauss_escalable(a,n)
    p = eye(n,n)
    l = zeros(n,n)
    if det(a) then
        for j=1:n-1
            [a,p,l]=permutacion_escalable(a,n,j,p,l)//Permuta las filas de las matriz
            for i=(j+1):n//Completa la matriz L y modifica las matriz A para llegar a una matriz triangular superior
                l(i,j) = a(i,j)/a(j,j)
                a=sumar_fila(a,i,j,l(i,j),n)
            end
        end
        for i=1:n//Agrega los 1 en la diagonal
           l(i,i)=1
        end
    else
        printf("Matriz singular\n");
    end
endfunction

function S = resolver_sin_pivote(A,n,vcons)
    //Resuelve un sistema de ecuaciones sin permutar filas
    S = zeros(n,1)
    [L,A] = gauss_sin_pivote(A,n)//Factoriza la matriz A=LU
    if verificar_matriz(L,n) then//Si no fue necesario intercambiar filas
        vcons = inv(L) * vcons'
        S(n)=vcons(n)/A(n,n)//Resuelvo por sustitucion hacia atras
        for i = 1:n-1
            s=0
            e=n-i
            for j= e+1:n
                s=s+A(e,j)*S(j)
            end
            S(e)=(vcons(e)-s)/A(e,e)
        end
    end
endfunction

function S = resolver_con_pivote(A,n,vcons)
    //Resuelve un sistema de ecuaciones permutando filas
    S = zeros(n,1)
    [P,L,A] = gauss_con_pivote(A,n)//Aplica la factorizacion PA=LU
    vcons = P * vcons'
    vcons = inv(L) * vcons
    S(n)=vcons(n)/A(n,n)
    for i = 1:n-1 //Resuelve haciendo sustitucion hacia atras
        s=0
        e=n-i
        for j= e+1:n
            s=s+A(e,j)*S(j)
        end
        S(e)=(vcons(e)-s)/A(e,e)
    end
endfunction

function S = resolver_escalable(A,n,vcons)
    //Resuelve un sistema de ecuaciones aplicacion la factorizacion descripta en el ejercicio 3
    S = zeros(n,1)
    [P,L,A] = gauss_escalable(A,n)//Factoriza la matriz
    vcons = P * vcons'
    vcons = inv(L) * vcons
    S(n)=vcons(n)/A(n,n)
    for i = 1:n-1//Resuelve aplicando sustitucion hacia atras
        s=0
        e=n-i
        for j= e+1:n
            s=s+A(e,j)*S(j)
        end
        S(e)=(vcons(e)-s)/A(e,e)
    end
endfunction
