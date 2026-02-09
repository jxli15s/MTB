
%%
e=1.602176634*10^(-19) % e=1.602176634*10^{-19} C
epsilon_0=8.854187817*10^(-12) % F/m
k_e=1/(4*pi*epsilon_0) %N*m2/C^2
Lm=7*10^(-9) % nm
unit=e^2*k_e*1/e % J to eV*A
%%
r=10^(-10);
d=4e-9;
epsilon=10;
z=10^6;
V=dual_gate_potential(r,d,epsilon,z)*unit %eV

%%
plot_dual_gate_potential(d,epsilon,unit)
%%
function plot_dual_gate_potential(d,epsilon,unit)
    % Constants
    e = 1.602176634e-19;        % Elementary charge (Coulombs)
    epsilon_0 = 8.854187817e-12; % Vacuum permittivity (F/m)
    k_e = 1 / (4 * pi * epsilon_0); % Coulomb constant (N·m²/C²)

    % Parameters
    % d = 1e-9;       % Distance to the gate (m)
    % epsilon = 3.9;  % Dielectric constant (e.g., SiO2)
    n_max = 1e6;    % Number of image charges to consider
    r = linspace(0.1e-10, 5e-10, 50); % Radial distance (m)

    % Calculate potential
    V = zeros(size(r)); % Initialize potential array
    for i = 1:length(r)
        V(i) = dual_gate_potential(r(i), d, epsilon, n_max)*unit;
    end

    % Plot the results
    figure;
    plot(r * 1e9, V, 'b-', 'LineWidth', 1.5); % Convert r to nm
    xlabel('Radial distance r (nm)', 'FontSize', 12);
    ylabel('Potential V (eV)', 'FontSize', 12);
    title('Dual-Gate Coulomb Potential vs Radial Distance', 'FontSize', 14);
    grid on;
end



function V = dual_gate_potential(r, d, epsilon, n_max)
    % Calculate dual-gate Coulomb potential using image charge method
    % r: radial distance
    % d: distance to the gates
    % q: charge
    % epsilon: dielectric constant
    % n_max: number of image charges to consider
    e=1.602176634*10^(-19);
    % epsilon_0=8.854187817*10^(-12);
    % k_e=1/(4*pi*epsilon_0);
    V = 0; % Initialize potential
    for n = -n_max:n_max
        V = V + ((-1)^n) / sqrt(r^2 + (2*n*d)^2);
    end
    % V = e^2/(4*pi*epsilon_0)/e/epsilon * V; % Final potential meV*m
    V=1/epsilon * V;
end