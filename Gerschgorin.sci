function [C,R] = Gerschgorin(A,n)
    C=%inf
    R=-%inf
    for i=1:n
        r=0
        d=0
        c=0
        for j=1:n
            if i==j then
                c=A(i,j)
            else
                r=r+abs(A(i,j))
            end
        end
        d=r+c
        if C>c then
            C=c
        end
        if R<d then
            R=d
        end
    end
    R=R-C
endfunction
