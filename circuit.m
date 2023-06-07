function CCTDesign=circuit(~)
% Define the component value ranges
Rin_range = [5000, 15000]; % Input impedance (ohm)
Rc_range = [500, 1500]; % Collector resistor (ohm)
RL_range = [5000, 15000]; % Load resistor (ohm)
Re_range = [500, 15000]; % Emitter resistor
fL_range = [10000, 50000]; % Lower cutoff frequency (Hz)
fH_range = [1000000, 20000000]; % Upper cutoff frequency (Hz)
Av_range = [50, 200]; % Voltage gain

num_tests = 100; % Number of random tests to perform

for test = 1:num_tests
    % Randomly select component values within the specified ranges
    Rin = randi(Rin_range);
    Rc = randi(Rc_range);
    RL = randi(RL_range);
    Re = randi(Re_range);
    fL = randi(fL_range);
    fH = randi(fH_range);
    Av = randi(Av_range);
    
    % Calculate the biasing circuit
    Vin = 0.1; % Input voltage (V)
    Vcc = 12; % Power supply voltage (V)
    Vb = Vcc*Rc/(Rc+RL); % Bias voltage (V)
    IcQ = Vcc/(2*Rc); % Quiescent collector current (A)
    Rb = (Vcc-Vb)/(Vin/IcQ - 1); % Base resistor (ohm)

    % Calculate the component values
    C1 = 1/(2*pi*fL*Rin); % Input coupling capacitor (F)
    C2 = 1/(2*pi*fH*RL); % Output coupling capacitor (F)
    Ce = (1 / (2 * pi * fL * Re)); % Emitter capacitance (F)
    Re = Vb/IcQ; % Emitter resistor (ohm)
    beta = 100; % DC current gain
    R1 = Rb/(beta+1); % Biasing resistor (ohm)
    R2 = Rb-R1; % Biasing resistor (ohm)

    % Create the circuit model
    gm = IcQ/(26*Re); % Transconductance (S)
    rpi = beta/gm; % Input resistance (ohm)
    Ro = RL; % Output resistance (ohm)
    Cpi = 2.6e-12; % Input capacitance (F)
    Cmu = 1.3e-12; % Miller capacitance (F)
    Hfe = beta/(1+beta); % Small signal current gain
    Av_calc = -gm*RL*Ro/(1+gm*Ro*(Cpi+Cmu)*Hfe); % Calculated voltage gain

    % Calculating gain difference
    desired_gain = 20*log10(Av); % desired gain in dB
    actual_gain = abs(20*log10(Av_calc));
    
    %Bandwidth
    Bw=fH-fL;
    
    % Calculate distortion
    Distortion = abs(100*(1-Hfe)/(2*Av_calc));
    
    % Calculate noise
    kTq = 1.38e-23*293; % Thermal voltage noise density at room temperature
    In = 2*beta*kTq/IcQ; % Input noise density (A/sqrt(Hz))
    En = In*Rin; % Input-referred noise (V/sqrt(Hz))
    Noise = abs(20*log10(En/(0.1*Av_calc)));
    
    CCTDesign=[Noise,Distortion];

end
%fprintf('Rin = %d ohm\n', Rin);
%fprintf('Rc = %d ohm\n', Rc);
%fprintf('Re = %d ohm\n', Re);
%fprintf('RL = %d ohm\n', RL);
%fprintf('fL = %d Hz\n', fL);
%fprintf('fH = %d Hz\n', fH);
end