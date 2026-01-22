function G= givens2x2(A)
    %Računanje Givensove rotacije za poništavanje donjeg lijevog elementa 2x2 matrice
    
    a = A(1, 1);
    b = A(2, 1);
    c = A(2, 1);
    d = A(2, 2);
    
    if abs(b) < 1e-16
        c = 1;
        s = 0;
        %posebni slučaj za simetrične matrice
    elseif abs(b - c) < 100*eps*max(1, norm(A, 'fro'));
        tau = (d - a)/(2*b);
        t = sign(tau) / (abs(tau)+ sqrt(1+tau^2));
        c = 1/sqrt(1 + t^2);
        s = t * c;
    else
        if abs(b) <= abs(a)
            z = b/a;
            c = 1/sqrt(1 + z^2);
            s = z/sqrt(1 + z^2);
        else
            z = a/b;
            c = abs(z)/sqrt(1 + z ^2);
            s = sign(z)/sqrt(1 + z^2);
        end
    end
    
    G = [c, s;-s, c];
    
end