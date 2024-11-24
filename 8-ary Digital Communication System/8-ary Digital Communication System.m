%% Plotting the SER curve
clear
N = 10^6; % Total number of symbols to transmit
Es_N0_dB = [0:18]; % Range of Es/N0 in dB
Es_N0 = 10.^(Es_N0_dB/10); % Convert to linear scale
Es = 6; % Choosing d=2 and using average energy formula 3/2*d^2 as derived 
N0 = Es./Es_N0; % Noise power spectral density

s1 = 2*randi([0,1],1,N) - 1; % phi(1) coordinate of signal vector contains binary PAM alphabets [-1,1]
s2 = 2*randi([0,3],1,N) - 3; % phi(2) coordinate of signal vector contains 4-PAM alphabets [-3,-1,1,3]

r1_hat = zeros(1,N); % initialising phi(1) coordinate of demodulated receive vector
r2_hat = zeros(1,N); % initialising phi(2) coordinate of demodulated receive vector

for i=1:length(Es_N0_dB)
    % Bit to symbol mapping
    n1 = sqrt(N0(i)/2).*randn(1,N); % phi(1) coordinate of white gaussian noise vector with variance = N0/2
    n2 = n1; % phi(2) coordinate of white gaussian noise vector with variance = N0/2
    r1 = s1 + n1; % phi(1) coordinate of receive vector
    r2 = s2 + n2; % phi(2) coordinate of receive vector

    % Demodulation/Symbol decisions
    r1_hat(r1 < 0) = -1; % assign -1 if phi(1) coordinate of receive vector is less than 0
    r1_hat(r1 >= 0) = 1; % assign 1 if phi(1) coordinate of receive vector is greater than 0 or equal to 0

    r2_hat(r2 < -2) = -3; % assign -3 if phi(2) coordinate of receive vector is less than -2
    r2_hat(r2 > 2) = 3; % assign 3 if phi(2) coordinate of receive vector is greater than 2
    r2_hat(r2 > -2 & r2 <= 0) = -1; % assign -1 if phi(2) coordinate of receive vector is between -2 and 0
    r2_hat(r2 > 0 & r2 <= 2) = 1; % assign 1 if phi(2) coordinate of receive vector is between 0 and 2
    
    % Symbol to bit mapping
    n_errors(i) = size(find((r1_hat - s1) | (r2_hat - s2)),2); % count the number of errors
end

simSER = n_errors/N; % find simulation SER by dividing the number of symbol errors by the total number of symbols

q_func = qfunc(sqrt(Es_N0/3)); % find qfunction value using similar erfc MATLAB command   also = (1/2)*erfc(sqrt(Es_N0/6))
simple_theoreticalSER = (5/2)*q_func; % simplified theoretical formula for SER
theoreticalSER = (5/2)*q_func - (3/2)*(q_func.^2); % actual theoretical formula for SER
figure
semilogy(Es_N0_dB,simSER,'go-'); % plot simulation SER against Es/N0 in dB
hold on;
semilogy(Es_N0_dB,simple_theoreticalSER,'mx-'); % plot simplified theoretical SER against Es/N0 in dB
hold on;
semilogy(Es_N0_dB,theoreticalSER,'rx-'); % plot actual theoretical SER against Es/N0 in dB
grid on
legend('Simulation', 'Simplified', 'Theoretical');
xlabel('Es/No (dB)')
ylabel('Symbol Error Rate (SER)')
title('Symbol error rate curve for 8-QAM')
%% Plotting the BER curve
clear
N = 10^6; % Total number of symbols to transmit
k = 3; % 3 bits per symbol
Eb_N0_dB = [0:12]; % Range of Eb/N0 in dB

Eb_N0 = 10.^(Eb_N0_dB/10); % Convert Eb/N0 to linear scale
Eb = 2; % Eb = Es/3 in this case
N0 = Eb./Eb_N0; % Noise power spectral density

symbols1 = 2*[0:1] - 1; % phi(1) coordinate of signal vector contains binary PAM alphabets [-1,1]
symbols2 = 2*[0:3] - 3; % phi(2) coordinate of signal vector contains 4-PAM alphabets [-3,-1,1,3]

r1_hat = zeros(1,N); % initialising phi(1) coordinate of demodulated receive vector
r2_hat = zeros(1,N); % initialising phi(2) coordinate of demodulated receive vector

for i=1:length(Eb_N0_dB)
    % Bit to symbol mapping
    blockBits = reshape((rand(1,N*k,1)<0.5),N,k); % 3-bit blocks of random bits
    sGrayCode1 = blockBits(:,1); % phi(1) coordinate gray code
    sBits2 = blockBits(:,[2:3]); % getting phi(2) coordinate bit values
    sDecimal2 = sum(sBits2.*[2 1],2); % converting binary to decimal format for gray code calculation
    sGrayCode2 = bitxor(sDecimal2,floor(sDecimal2/2)); % phi(2) coordinate gray codes
    
    s1 = symbols1(sGrayCode1+1); % phi(1) coordinate of signal vector using gray code mapping
    s2 = symbols2(sGrayCode2+1); % phi(2) coordinate of signal vector using gray code mapping

    n1 = sqrt(N0(i)/2).*randn(1,N); % phi(1) coordinate of white gaussian noise vector with variance = N0/2
    n2 = n1; % phi(2) coordinate of white gaussian noise vector with variance = N0/2
    r1 = s1 + n1; % phi(1) coordinate of receive vector
    r2 = s2 + n2; % phi(2) coordinate of receive vector
    
    % Demodulation/Symbol decisions
    r1_hat(r1 < 0) = -1; % assign -1 if phi(1) coordinate of receive vector is less than 0
    r1_hat(r1 >= 0) = 1; % assign 1 if phi(1) coordinate of receive vector is greater than 0 or equal to 0

    r2_hat(r2 < -2) = -3; % assign -3 if phi(2) coordinate of receive vector is less than -2
    r2_hat(r2 > 2) = 3; % assign 3 if phi(2) coordinate of receive vector is greater than 2
    r2_hat(r2 > -2 & r2 <= 0) = -1; % assign -1 if phi(2) coordinate of receive vector is between -2 and 0
    r2_hat(r2 > 0 & r2 <= 2) = 1; % assign 1 if phi(2) coordinate of receive vector is between 0 and 2
    
    % Symbol to bit mapping
    GrayCodeLUT2 = [0 0; 0 1; 1 1; 1 0]; % Gray code look-up table for converting from symbols to bits
    rGrayCode1 = floor((r1_hat+2)/2); % Converts symbols [-1,1] to 1-bit gray code [0,1]
    rGrayCode2 = GrayCodeLUT2(floor((r2_hat+4)/2+1),:); % Converts symbols [-3,-1,1,3] to 2-bit gray code
    n_errors(i) = size(find([rGrayCode1' rGrayCode2] - blockBits),1); % count the number of bit errors
end

simBER = n_errors/(N*k); % find simulation BER by dividing the number of bit errors by the total number of bits transmitted

q_func = qfunc(sqrt(Eb_N0)); % find qfunction value using similar erfc MATLAB command   also = (1/2)*erfc(sqrt(Eb_N0/2))
theoreticalBER = (5/6)*q_func - (1/2)*(q_func.^2); % simplified theoretical formula for BER
figure
semilogy(Eb_N0_dB,simBER,'go-'); % plot simulation BER against Eb/N0 in dB
hold on;
semilogy(Eb_N0_dB,theoreticalBER,'mx-'); % plot theoretical BER against Eb/N0 in dB
grid on
legend('Simulation', 'Theoretical');
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate (BER)')
title('Bit error rate curve for 8-QAM')