function e = eigs_from_schur(T)
    %Čitanje svojstvenih vrijednosti iz realne Schurove forme
    n = length(T);
    
    e = zeros(n, 1);
    
    i = 1;
    while(i <= n)
    
        if i == n %ako smo došli u posljednji element on mora biti realna svojstvena vrijednosti
            e(i) = T(i, i);
            i = i+1;
        elseif (T(i+1, i) == 0) %ako je poddijagonalni element = 0, T(i+1, i) je svojstvena vrijednost
            e(i) = T(i, i);
            i  = i+1;
            else %ako je poddijagonalni element ~= 0, svojestvene vrijednosti dobijemo iz 2x2 bloka
            tr = T(i, i) + T(i+1, i+1);
            det = T(i, i)*T(i+1, i+1) - T(i, i+1)*T(i+1, i);
            disc = tr^2 - 4*det;
            
            if tr >= 0
                e(i) = (tr + sqrt(disc))/2;
            else
                e(i) = (tr - sqrt(disc))/2;
            end
            e(i+1) = det/e(i);
            
            i = i + 2;
        end
    
    end
    
 end