n = 10;

A = randn(n);

[T, Q, it, f] = real_schur(A);

% Crtamo realnu Schurovu formu. Crne točke označavaju elemente jednake 0.
figure(1);
imagesc(log10(abs(T))); colorbar; axis square; hold on;
[xind, yind] = find(T == 0);
mark_zero = scatter(yind, xind, 100, 'k', 'filled'); hold off;
title("log10(abs(T))");
legend(mark_zero, 'T == 0', 'FontSize', 18);

if f == 1
    fprintf("Algoritam je konvergirao u %i koraka.\n\n", it);
else
    fprintf("Algoritam nije konvergirao.\n\n");
end

fprintf("||T - Q'*A*Q||: %.5e\n", norm(T - Q'*A*Q, 'fro'));
fprintf("||I - Q'*Q||: %.5e\n\n", norm(eye(n) - Q'*Q, 'fro'));

original_eigs = sort(eig(A));
schur_eigs = sort(eigs_from_schur(T));

fprintf("Suma razlika svojstvenih vrijednosti: %.5e\n\n", norm(original_eigs - schur_eigs, 1));

fprintf("Sortirane svojstvene vrijednosti:\n");
fprintf("eig(A)               eig_from_schur(T)\n");
for k = 1:length(original_eigs)
    fprintf("|% .5f %+ .5fi  | % .5f %+ .5fi|\n", ...
        real(original_eigs(k)), imag(original_eigs(k)), ...
        real(schur_eigs(k)), imag(schur_eigs(k)));
end