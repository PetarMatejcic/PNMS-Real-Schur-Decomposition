function plot_schur_state(H, p, q, pause_s, small_tol)
%PLOT_SCHUR_STATE Visualize Hessenberg matrix during real Schur iteration
%
%   plot_schur_state(H, p, q, pause_s, small_tol)
%
%   H         - current Hessenberg matrix
%   p, q      - active block is rows/cols p+1:q
%   pause_s   - pause duration in seconds (0 for no pause)
%   small_tol - threshold for "numerically small" (e.g. 1e-10)

    if nargin < 5
        small_tol = 1e-10;
    end
    if nargin < 4
        pause_s = 0;
    end

    n = size(H,1);

    % ---- magnitude plot (log10) ----
    clf;
    mag = abs(H);
    mag(mag == 0) = NaN;    % true zeros appear as white
    imagesc(log10(mag));
    axis equal;
    colormap(jet);
    colorbar;
    hold on;

    % ---- mark numerically small but nonzero entries ----
    [i_small, j_small] = find(abs(H) > 0 & abs(H) < small_tol);
    plot(j_small, i_small, 'k.', 'MarkerSize', 6)

    % ---- draw active block box ----
    if q >= p+1
        x = [p+0.5, q+0.5, q+0.5, p+0.5, p+0.5];
        y = [p+0.5, p+0.5, q+0.5, q+0.5, p+0.5];
        plot(x, y, 'r-', 'LineWidth', 2)
    end

    % ---- draw Hessenberg subdiagonal for reference ----
    plot(1:n-1, 2:n, 'w:', 'LineWidth', 1)

    hold off

    title(sprintf('log_{10}(|H|), active block = [%d:%d]', p+1, q))
    xlabel('column')
    ylabel('row')

    drawnow

    if pause_s > 0
        pause(pause_s)
    end
end
