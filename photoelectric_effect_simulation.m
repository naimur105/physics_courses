% Photoelectric Effect Simulation
% Author: [Naimur Rahman, ChatGPT]

clc;
clear;
close all;

% Constants
h = 6.626e-34; % Planck's constant (J·s)
c = 3e8;       % Speed of light (m/s)
e = 1.602e-19; % Electron charge (C)
m_e = 9.109383e-31

% Input Parameters
fprintf('Photoelectric Effect Simulation\n');
lambda = input('Enter wavelength (nm): ') * 1e-9; % Convert nm to meters
intensity = input('Enter intensity (W/m^2): ');  % Light intensity
phi = input('Enter work function (eV): ') * e;   % Convert eV to Joules
V_range = input('Enter voltage range as [V_min V_max]: '); % Voltage range for V-I plot

% Derived Parameters
frequency = c / lambda; % Frequency of light
energy_per_photon = h * frequency; % Energy of a single photon
fprintf('Photon energy: %.3f eV\n', energy_per_photon / e);

% Check if photoemission occurs
if energy_per_photon < phi
    error('Photon energy is less than the work function. No photoemission occurs.');
end

% Compute photoelectric current
V = linspace(V_range(1), V_range(2), 500); % Voltage range
KE_max = energy_per_photon - phi; % Maximum kinetic energy of emitted electrons (J)
V_stop = KE_max / e; % Stopping potential (V)

% Electron emission rate (current proportional to intensity)
electron_rate = intensity / energy_per_photon; % Number of photons per second per unit area
current = zeros(size(V)); % Initialize current array
for i = 1:length(V)
    if V(i) <= V_stop
        KE = e * (KE_max - e * V(i)); % Kinetic energy decreases with voltage
        current(i) = max(0, e * electron_rate * sqrt(KE / (2 * m_e))); % Current proportional to KE
    else
        current(i) = 0; % No current beyond stopping potential
    end
end

% Plot V-I Curve
figure;
plot(V, current, 'LineWidth', 2);
xlabel('Voltage (V)');
ylabel('Current (A)');
title('Photoelectric Effect: V-I Curve');
grid on;

% V-I Plot for Different Intensities
intensities = [intensity, intensity * 2, intensity * 3]; % Different intensities
figure;
hold on;
for k = 1:length(intensities)
    electron_rate = intensities(k) / energy_per_photon; % Update for new intensity
    current = zeros(size(V)); % Reset current
    for i = 1:length(V)
        if V(i) <= V_stop
            KE = e * (KE_max - e * V(i));
            current(i) = max(0, e * electron_rate * sqrt(KE / (2 * m_e)));
        else
            current(i) = 0; % No current beyond stopping potential
        end
    end
    plot(V, current, 'LineWidth', 2, 'DisplayName', sprintf('Intensity = %.1f W/m^2', intensities(k)));
end
xlabel('Voltage (V)');
ylabel('Current (A)');
title('Photoelectric Effect: V-I Curve with Different Intensities');
legend;
grid on;

% V-I Plot for Different Wavelengths (Frequencies)
wavelengths = [lambda, lambda * 0.9, lambda * 0.8]; % Shorter wavelengths (higher frequencies)
figure;
hold on;
for k = 1:length(wavelengths)
    frequency = c / wavelengths(k);
    energy_per_photon = h * frequency;
    KE_max = energy_per_photon - phi; % Update KE_max
    if KE_max > 0
        V_stop = KE_max / e; % Update stopping potential
        electron_rate = intensity / energy_per_photon; % Update for new wavelength
        current = zeros(size(V)); % Reset current
        for i = 1:length(V)
            if V(i) <= V_stop
                KE = e * (KE_max - e * V(i));
                current(i) = max(0, e * electron_rate * sqrt(KE / (2 * m_e)));
            else
                current(i) = 0; % No current beyond stopping potential
            end
        end
        plot(V, current, 'LineWidth', 2, 'DisplayName', sprintf('Wavelength = %.1f nm', wavelengths(k) * 1e9));
    else
        fprintf('No photoemission for wavelength = %.1f nm (photon energy < work function).\n', wavelengths(k) * 1e9);
    end
end
xlabel('Voltage (V)');
ylabel('Current (A)');
title('Photoelectric Effect: V-I Curve with Different Wavelengths');
legend;
grid on;

% V vs. frequency plot
frequencies = linspace(c / (lambda * 1.5), c / (lambda * 0.5), 500); % Frequency range for the plot
cutoff_frequency = phi / h; % Cutoff frequency
V_f = zeros(size(frequencies)); % Initialize voltage array

for i = 1:length(frequencies)
    if frequencies(i) >= cutoff_frequency
        energy_per_photon = h * frequencies(i); % Photon energy
        KE_max = energy_per_photon - phi; % Maximum kinetic energy
        V_f(i) = KE_max / e; % Voltage proportional to KE
    else
        V_f(i) = 0; % No photoemission below cutoff frequency
    end
end

% Plot V vs. frequency
figure;
plot(frequencies, V_f, 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Voltage (V)');
title('Photoelectric Effect: Voltage vs. Frequency');
grid on;

