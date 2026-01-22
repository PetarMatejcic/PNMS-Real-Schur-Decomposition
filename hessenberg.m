function [A, U] = hessenberg(A)
    %Računanje Hessenbergove forme
    
    n = length(A);
    U = eye(n);
    
    %Uzastopno poništavnje elementa ispod poddijagonale Householderovim reflektorima
    for k = 1:n-2
        
        [v, b] = house(A(k+1:n,k));
        P = eye(n-k) - b*v*v';
        A(k+1:n, k:n) = P*A(k+1:n, k:n);
        A(1:n, k+1:n) = A(1:n, k+1:n)*P;
        U(1:n, k+1:n) = U(1:n, k+1:n)*P;
        
    end
    
    %Postavljenje nula ispod dijagonale za ljepši ispis
    A = triu(A, -1);
    
end