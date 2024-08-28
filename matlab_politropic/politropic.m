
rho_cp_T1 = readmatrix("averages_heat_NPT1.dat", "Range", "F1:F34");
cp_T1 = readmatrix("averages_heat_NPT1.dat", "Range", "H1:H34");
rho_cp_T2 = readmatrix("averages_heat_NPT2.dat", "Range", "F1:F33");
cp_T2 = readmatrix("averages_heat_NPT2.dat", "Range", "H1:H33");


rho_cv_T1 = readmatrix("averages_heat_NVT1.dat", "Range", "F1:F34");
cv_T1 = readmatrix("averages_heat_NVT1.dat", "Range", "G1:G34");
rho_cv_T2 = readmatrix("averages_heat_NVT2.dat", "Range", "F1:F33");
cv_T2 = readmatrix("averages_heat_NVT2.dat", "Range", "G1:G33");


rhos = linspace(0.002, 1.3, 3000);
gammas1 = [];
for r = rhos
    cp = linearSampling(rho_cp_T1, cp_T1, [r]);
    cv = linearSampling(rho_cv_T1, cv_T1, [r]);
    gammas1(length(gammas1) + 1) = cp(1) / cv(1);
end
gammas1 = gammas1';

rhos = linspace(0.002, 1.3, 3000);
gammas2 = [];
for r = rhos
    cp = linearSampling(rho_cp_T2, cp_T2, [r]);
    cv = linearSampling(rho_cv_T2, cv_T2, [r]);
    gammas2(length(gammas2) + 1) = cp(1) / cv(1);
end
gammas2 = gammas2';

plot(rhos, gammas1, "DisplayName", "T = 1.0", 'LineWidth', 1.25, "Color", "blue");
hold on
plot(rhos, gammas2, "DisplayName", "T = 2.0", 'LineWidth', 1.25, "Color", "red");
grid on
grid minor
legend('show');
xlabel("Densità")
ylabel("γ = c_P/c_V")
xlim([0.0, 1.3])
title("Indice politropico")
% set(gca, 'YScale', 'log');

exportFigure("politropic.png");