n = 10;

A = randn(n);

[T, Q] = real_schur(A);

figure(1);
imagesc(log10(abs(T))); colorbar; axis square;

fprintf("||T - Q'*A*Q||: %.5e\n", norm(T - Q'*A*Q, 'fro'));
fprintf("||I - Q'*Q||: %.5e\n\n", norm(eye(n) - Q'*Q, 'fro'));

original_eigs = sort(eig(A));
schur_eigs = sort(eigs_from_schur(T));

fprintf("Sortirane svojstvene vrijednosti:\n");
fprintf("eig(A)               eig_from_schur(T)\n");
for k = 1:length(original_eigs)
    fprintf("|% .5f %+ .5fi  | % .5f %+ .5fi|\n", ...
        real(original_eigs(k)), imag(original_eigs(k)), ...
        real(schur_eigs(k)), imag(schur_eigs(k)));
end