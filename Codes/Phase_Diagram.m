% Name: Mohammad Reza Hadi
% E-mail: mohammadreza.hadi.2002@gmail.com

clc; clear; close all;

%% ===================== CONSTANTS =====================
T_m_A = 1135;        % K
T_m_B = 1685;        % K

delta_H_A = 19.2e3;  % J/mol
delta_H_B = 22.8e3;  % J/mol

Omega_S = 13.9e3;    % J/mol
Omega_L = 9.1e3;     % J/mol

R = 8.3144;          % J/mol.K

%% ===================== GRIDS =====================
T = 700:1:1800;     % temperature grid
eps = 1e-6;
X_B = eps:0.001:1-eps;
X_A = 1 - X_B;

%% ===================== STORAGE =====================
X_S1 = nan(size(T));   % solid low
X_S2 = nan(size(T));   % solid high
X_L1 = nan(size(T));   % liquid low
X_L2 = nan(size(T));   % liquid high

%% ===================== SOLVER OPTIONS =====================
options = optimoptions('fsolve','Display','off','FunctionTolerance',1e-12,'StepTolerance',1e-12);
X_tol = 1e-6;   % minimum ΔX to consider two-phase coexistence

%% ===================== MAIN LOOP =====================

for k = 1:length(T)

    T0 = T(k);

    % Gibbs energy of melting for A and B
    dGm_A = delta_H_A * ((T_m_A - T0) / T_m_A);
    dGm_B = delta_H_B * ((T_m_B - T0) / T_m_B);

    % Gibbs energies curves for solid and liquid
    G_S = -X_A .* dGm_A + R*T0*(X_A.*log(X_A) + X_B.*log(X_B)) + Omega_S .* X_A .* X_B;
    G_L =  X_B .* dGm_B + R*T0*(X_A.*log(X_A) + X_B.*log(X_B)) + Omega_L .* X_A .* X_B;

    % Numerical derivatives
    dG_S = gradient(G_S, X_B);
    dG_L = gradient(G_L, X_B);

    % Interpolants
    GS  = @(x) interp1(X_B, G_S, x, 'pchip');
    GL  = @(x) interp1(X_B, G_L, x, 'pchip');
    dGS = @(x) interp1(X_B, dG_S, x, 'pchip');
    dGL = @(x) interp1(X_B, dG_L, x, 'pchip');
    
    % Solid minima
    [~, locs_S] = findpeaks(-G_S);
    % Liquid minima
    [~, locs_L] = findpeaks(-G_L);

    min_XS = X_B(locs_S);
    min_XL = X_B(locs_L);

    found = false;

    try
        % ----------------- SOLID-LIQUID TANGENT -----------------
        % Only attempt if Gibbs curves cross
        if any(G_L - G_S > 0) && any(G_L - G_S < 0)
            
            % Better x0: intersection estimate
            [~, idx_cross] = min(abs(G_L - G_S));
            x0_guess = X_B(idx_cross);
            if (T0 >= T_m_A) && (T0 < 1151)
                for i = 1:length(min_XS)
                    for j = 1:length(min_XL)
                        
                        x0 = [0.01 0.02];
    
                        F = @(x) [dGS(x(1))-dGL(x(2)); (GL(x(2))-GS(x(1)))/(x(2)-x(1))-dGS(x(1))];
                        xeq = fsolve(F, x0, options);
                        
                        if all(xeq>0 & xeq<1) && abs(xeq(2)-xeq(1))>X_tol
                            X_S2(k) = xeq(1);
                            X_L1(k) = xeq(2);
                            found = true;
                            break
                        end
                    end
                    if found, break; end
                end
            
            else
                for i = 1:length(min_XS)
                    for j = 1:length(min_XL)
                        % Start near the intersection and minima
                        x0 = [min_XS(i)*(1-0.2)+x0_guess*0.2, min_XL(j)*(1-0.2)+x0_guess*0.2];
    
                        F = @(x) [dGS(x(1))-dGL(x(2)); (GL(x(2))-GS(x(1)))/(x(2)-x(1))-dGS(x(1))];
                        xeq = fsolve(F, x0, options);
                        
                        if all(xeq>0 & xeq<1) && abs(xeq(2)-xeq(1))>X_tol
                            X_S2(k) = xeq(1);
                            X_L1(k) = xeq(2);
                            found = true;
                            break
                        end
                    end
                    if found, break; end
                end
            end

        % ----------------- SOLID-SOLID TANGENT -----------------
        elseif any(diff(sign(gradient(G_S)))~=0)
            F = @(x) [dGS(x(1))-dGS(x(2)); (GS(x(2))-GS(x(1)))/(x(2)-x(1))-dGS(x(1))];
            x0 = [0.2 0.8];
            xeq = fsolve(F, x0, options);
            if abs(xeq(2)-xeq(1)) > X_tol
                X_S1(k) = xeq(1);
                X_S2(k) = xeq(2);
            end

        % ----------------- LIQUID-LIQUID TANGENT -----------------
        elseif any(diff(sign(gradient(G_L)))~=0)
            F = @(x) [dGL(x(1))-dGL(x(2)); (GL(x(2))-GL(x(1)))/(x(2)-x(1))-dGL(x(1))];
            x0 = [0.2 0.8];
            xeq = fsolve(F, x0, options);
            if abs(xeq(2)-xeq(1)) > X_tol
                X_L1(k) = xeq(1);
                X_L2(k) = xeq(2);
            end
        end

    catch
        % Skip if solver fails
    end
end

%% ===================== PHASE DIAGRAM =====================
figure
hold on
plot(X_S1, T, 'r', 'LineWidth', 2)
plot(X_L1, T, 'b', 'LineWidth', 2)
plot(X_S2, T, 'r', 'LineWidth', 2)
plot(X_L2, T, 'b', 'LineWidth', 2)
xlabel('X_B')
ylabel('Temperature (K)')
title('Binary Phase Diagram', ...
    'Interpreter','latex',...
    'FontSize',16)
text(0.42, 750, '$\alpha_1 + \alpha_2$', ...
     'Interpreter','latex', ...
     'FontSize',16)
text(0.42, 1300, 'S + L', ...
     'Interpreter','latex', ...
     'FontSize',16)
text(0.49, 1000, 'S', ...
     'Interpreter','latex', ...
     'FontSize',16)
text(0.49, 1550, 'L', ...
     'Interpreter','latex', ...
     'FontSize',16)
fontname("Times New Roman")
legend('Solidus','Liquidus','Location','northwest')
grid on
hold off

%% ===================== SAMPLE COMMON-TANGENT PLOT =====================
T_samples = [700,750,800,1200,1400,1600];   % six temperatures

figure('WindowState','maximized')
t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
title(t,'Common-Tangent Points: ($\Delta G^M_L$ and $\Delta G^M_S$) vs $X_B$', ...
      'Interpreter','latex', ...
      'FontSize',22)

for k = 1:length(T_samples)

    T_plot = T_samples(k);
    idx = find(T == T_plot, 1);

    nexttile
    hold on
    grid on
    fontname("Times New Roman")

    if isempty(idx)
        title(['T = ' num2str(T_plot) ' K (no data)'])
        continue
    end

    T0 = T_plot;
    dGm_A = delta_H_A * ((T_m_A - T0)/T_m_A);
    dGm_B = delta_H_B * ((T_m_B - T0)/T_m_B);

    G_S = -X_A .* dGm_A + R*T0*(X_A.*log(X_A)+X_B.*log(X_B)) + Omega_S.*X_A.*X_B;
    G_L =  X_B .* dGm_B + R*T0*(X_A.*log(X_A)+X_B.*log(X_B)) + Omega_L.*X_A.*X_B;

    dG_S = gradient(G_S, X_B);
    dG_L = gradient(G_L, X_B);

    GS  = @(x) interp1(X_B, G_S, x, 'pchip');
    GL  = @(x) interp1(X_B, G_L, x, 'pchip');
    dGS = @(x) interp1(X_B, dG_S, x, 'pchip');
    dGL = @(x) interp1(X_B, dG_L, x, 'pchip');

    % ---- Plot Gibbs energies ----
    plot(X_B, G_S, 'r', 'LineWidth', 2)
    plot(X_B, G_L, 'b', 'LineWidth', 2)

    % ---- Plot common tangent ----
    if ~isnan(X_S1(idx)) && ~isnan(X_S2(idx))
        x1 = X_S1(idx); x2 = X_S2(idx);
        m = dGS(x1); b = GS(x1) - m*x1;
        plot([x1 x2], m*[x1 x2] + b, 'k--', 'LineWidth', 2)
        plot(x1, GS(x1), 'ro', 'MarkerSize', 7, 'LineWidth', 2)
        plot(x2, GS(x2), 'ro', 'MarkerSize', 7, 'LineWidth', 2)
        
        drawnow
        yl = ylim;
        plot([x1 x1], [yl(1) GS(x1)], 'k:', 'LineWidth',1.5)
        plot([x2 x2], [yl(1) GS(x2)], 'k:', 'LineWidth',1.5)

    elseif ~isnan(X_S2(idx)) && ~isnan(X_L1(idx))
        x1 = X_S2(idx); x2 = X_L1(idx);
        m = dGS(x1); b = GS(x1) - m*x1;
        plot([x1 x2], m*[x1 x2] + b, 'k--', 'LineWidth', 2)
        plot(x1, GS(x1), 'ro', 'MarkerSize', 7, 'LineWidth', 2)
        plot(x2, GL(x2), 'bo', 'MarkerSize', 7, 'LineWidth', 2)
        
        drawnow
        yl = ylim;
        plot([x1 x1], [yl(1) GS(x1)], 'k:', 'LineWidth',1.5)
        plot([x2 x2], [yl(1) GL(x2)], 'k:', 'LineWidth',1.5)

    elseif ~isnan(X_L1(idx)) && ~isnan(X_L2(idx))
        x1 = X_L1(idx); x2 = X_L2(idx);
        m = dGL(x1); b = GL(x1) - m*x1;
        plot([x1 x2], m*[x1 x2] + b, 'k--', 'LineWidth', 2)
        plot(x1, GL(x1), 'bo', 'MarkerSize', 7, 'LineWidth', 2)
        plot(x2, GL(x2), 'bo', 'MarkerSize', 7, 'LineWidth', 2)
        
        drawnow
        yl = ylim;
        plot([x1 x1], [yl(1) GL(x1)], 'k:', 'LineWidth',1.5)
        plot([x2 x2], [yl(1) GL(x2)], 'k:', 'LineWidth',1.5)
    end

    xlabel('X_B',...
        'FontSize',20)
    ylabel('\DeltaG (J/mol)',...
        'FontSize',20)
    title(['T = ' num2str(T_plot) ' K'],...
        'FontSize',16)
    lgd = legend('\DeltaG^M_S','\DeltaG^M_L','Common Tangent','Location','north','fontsize',12);
    lgd.NumColumns = 3;
end
