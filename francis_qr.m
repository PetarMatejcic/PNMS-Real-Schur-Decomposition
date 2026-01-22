function [H, Z] = francis_qr(H)
    %Implicitni QR algoritam s pomacima
    
    n = length(H);
    m = n - 1;
    
    %Računanje koefcijenta karakterističnog polinoma donjeg desnog 2x2 bloka
    s = H(m, m) + H(n, n);
    t = H(m, m)*H(n, n) - H(m, n)*H(n, m);

    %Racunanje prvog stupaca matrice (H^2 - sH +tI)
    %Za Hessenbergovu matricu on ima samo prva tri elementa ~= 0 pa samo njih racunamo
    x = H(1, 1)^2 + H(1, 2)*H(2, 1) - s*H(1, 1) + t;
    y = H(2, 1)*(H(1, 1) + H(2, 2) - s);
    z = H(2, 1)*H(3, 2);
    
    %Inicijalizacija matrice transformacija
    Z = eye(n);
    
    for k = 0:n-3
    %Računanje Householderovog reflektora koji poništava prvi stupac (y i z)
        [v, b] = house([x;y;z]);
        P = eye(3) - b*v*v';
        
        %Primijena P' (odnosto P) na H s lijeva (odnosno desna)
        %Time nastaje "bulge", ali se očuva sličnost
        q = max(1, k);
        H(k+1:k+3, q:n) =P'*H(k+1:k+3,q:n);
        r = min(k+4, n);
        H(1:r, k+1:k+3) = H(1:r, k+1:k+3)*P;
        
        %Definiranje sljedećeg stupca koji će se poništiti
        x = H(k+2, k+1);
        y = H(k+3, k+1);
        if k < n-3
            z = H(k+4, k+1);
        end
        
        %Spremanje transformacije
        Z(1:r, k+1:k+3) = Z(1:r, k+1:k+3)*P;
        
        %U sljedećoj iteraciji poništavanjem sljedećeg stupca pomičemo "bulge" dolje desno.
        %Ovaj postupak ponavljamo dok jedini dio "bulge"-a ne bude element H(n, n-1).
        
    end
    
    %Poništavanje elementa H(n, n-1)
    [v, b] = house([x;y]);
    P = eye(2) - b*v*v';
    H(n-1:n, n-2:n) = P*H(n-1:n, n-2:n);
    H(1:n, n-1:n) = H(1:n, n-1:n)*P;
    Z(1:n, n-1:n) = Z(1:n, n-1:n)*P;
    
end