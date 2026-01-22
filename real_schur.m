function [H, Q, it, converge_flag] = real_schur(A, maxit = 1000, tol = 100*eps)
    
    % Posebni slučajevi (nekvadratna matrica, prazna matrica, skalar)
    [m, n] = size(A);
    if m != n
        error('Matrica nije kvadratna!');
        return;
    elseif n == 0
        error('Matrica je prazna');
        return;
    elseif n == 1
        H = A; Q = 1;
        return;
    end

    % Računanje Hessenbergove forme
    [H, Q] = hessenberg(A);
    
    p = 0;
    q = n;
    
    % Norma poddijagonale (za detekciju stagnacije)
    % Maksimalni broj iteracija na bloku i brojač iteracija unutar bloka
    block_maxit = 30 * max(10, q-p);
    block_it = 0;
    
    %Brojač ukupnog broja iteracija
    it = 0;
    
    while (q > 1 && it < maxit)
        it = it + 1;
        
        % Postavljanje relativno malih podijagonalnih elemenata na 0
        % Relativno mali u smislu da su mali u odnosu na susjedne dijagonale elemente i zadanu toleranciju.
        for i = 2:n
            if abs(H(i, i-1)) <= tol*(abs(H(i-1, i-1)) + abs(H(i, i)))
                H(i, i-1) = 0;
            end
        end
        
        % Određivanje donje granice aktivnog bloka tako što krećemo od dna trenutnog aktivnog bloka
        % i smanjujemo q dok ne dodjemo do podijagonalnog elementa ~= 0.
        while q>2 && H(q, q-1) == 0
            q = q - 1;
            block_it = 0;
            block_maxit = max(10, q - p);
        end
        
        % Određivanje gornje granice aktivnog bloka tako što krećemo od q-1 
        % i povećavamo p dok ne dođemo do poddijegonalnog elementa = 0
        p = q-1;
        while p >= 1 && H(p+1,p) ~= 0
            p = p-1;
            block_it = 0;
            block_maxit = max(10, q - p);
        end
        % Ako je aktivni blok duljine jedan, nalazimo se u 1x1 dijagonalnom bloku Schurove forma,
        % pa ga preskačemo i vraćamo se na početak while petlje
        if (q - p == 1)
            q = q - 1;
            it = 0;
            block_it = 0;
            block_maxit = max(10, q - p);
            continue;
        end
        
        subdiag_norm = norm(diag(H(p+1:q, p+1:q), -1), 1);
        
        % Ovisno o tome je li dimenzija aktivnog bloka =2 ili >2, drugačije postupamo
        if q - p == 2
            % Ako je aktivni blok dimenzije 2, provjeravamo ima li on kompleksne svojestvene vrijednosti
            % Ako ima prelazimo na sljedeći blok (trenutni aktivni blok je 2x2 blok Schurove forme)
            if (H(p+1, p+1) - H(q, q))^2 + 4*H(p+1, q)*H(q, p+1) < 0
                q = q-2;
            % U suprotnom, poništavamo donji ljevi element Givensovom rotacijom kako bi dobili gornje trokutastu matricu
            else
                G = givens2x2(H(p+1:q, p+1:q));
                H(:, p+1:q) = H(:, p+1:q) * G;
                H(p+1:q, :) = G' * H(p+1:q, :);
                Q(:, p+1:q) = Q(:, p+1:q) * G;
            end
            
        elseif q - p > 2
            % Izvođenje francis_qr na aktivnom bloku
            [H22, Z] = francis_qr(H(p+1:q, p+1:q));
            
            % Transformiranje matrice H i Q
            Q(:, p+1:q) = Q(:, p+1:q)*Z;
            H(:, p+1:q) = H(:, p+1:q)*Z;
            H(p+1:q, :) = Z'*H(p+1:q, :);
            H(p+1:q, p+1:q) = H22;
            
            % Provjera je li se norma poddijagonale značajno promijenila od prethodne
            new_subdiag_norm = norm(diag(H(p+1:q, p+1:q), -1), 1);
            if new_subdiag_norm/subdiag_norm >= 0.99
                block_it = block_it + 1;
            else
                block_it = 0;
            end
            
        end
        
        % Ako se norma poddijagonale nije značajno promijenila block_maxit iteracija,  izvršavamo izvanredan pomak
        % On je veći od običnih pomaka kako bi matricu izbacio iz stagnacije algoritma
        if block_it > block_maxit
            [H22, Z] = francis_qr(H(p+1:q, p+1:q), 0);
            Q(:, p+1:q) = Q(:, p+1:q)*Z;
            H(:, p+1:q) = H(:, p+1:q)*Z;
            H(p+1:q, :) = Z'*H(p+1:q, :);
            H(p+1:q, p+1:q) = H22;
        end
        
    end
    
    % Zadnje postavljanje nula
    for i = 2:n
        if abs(H(i, i-1)) <= tol*(abs(H(i-1, i-1)) + abs(H(i, i)))
            H(i, i-1) = 0;
        end
    end
    % Elemente ispod poddijagonale postavljamo na 0 (za ljepši ispis, a stvarni elementi su ~1e-16)
    H = triu(H, -1);
    
    % Oznaka je li algoritam konvergirao
    if it == maxit
        converge_flag = 0;
    else
        converge_flag = 1;
    end
    
end