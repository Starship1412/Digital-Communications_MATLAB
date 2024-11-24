% Defining variables
clear
d = 1; % chosen value of d, distance between the progin and symbol
k = 4; % message bits
n = 7; % number of bits of coded message
% N = 4*10^7; % length of data vector, total symbols tranmitted is n*N, data bits is k*N
N = 1e7; % For quick simulation
R = k/n; % Code Rate

%EbN0_dB = [0:0.1:10]; % Range of Eb/N0
Eb_N0_dB = [0:1:10]; % For quick simulation
Eb_N0 = 10.^(Eb_N0_dB/10); % in linear
Ec_N0 = Eb_N0*R;
Ec = d^2; % square of the euclidean distance from origin
N0 = Ec./Ec_N0;

P = [1 1 0; 0 1 1; 1 1 1; 1 0 1]; 
G = [P eye(k)]; % Generator matrix
H = [eye(n-k) transpose(P)]; % Parity Check Matrix
H_transpose = transpose(H);
error_index = [0; 3; 2; 5; 1; 7 ;4; 6]; % index where bit errors occur

u = rand(N,k) <0.5; % u, the message bits
% Encode using Hamming encoder
cBlocks = mod(u*G,2); % Coded message, using modulo-2 arithmetic

% BPSK Modulator
% Blocks of Data
signalBlock = zeros(size(cBlocks));
signalBlock(cBlocks==0) = -d;
signalBlock(cBlocks==1) = d;

% Upto here should be same for all Eb/N0, next addition of noise
% Go through for loop for each value of Eb/N0
n_errors = zeros(length(Eb_N0_dB),1); % pre-allocation for speed
sigma = sqrt(N0/2); %AWGN with variance = N0/2
rBitBlock = zeros(size(signalBlock)); % Allocating memory
for i = 1:length(Eb_N0_dB)
    % Add AWGN with variance = N0/2
    noiseBlock = randn(size(signalBlock)).*sigma(i);
    % Add noise to the transmit signal to create recieved signal
    recievedBlocks = signalBlock+noiseBlock;
    
    % Demodulatar (Hard Decisions)
    % recievedBlocks, contains -d,d + some noise
    % then if recievedBlocks is less than 0, then 0 is transmitted, else 1
    rBitBlock(recievedBlocks<0) = 0;
    rBitBlock(recievedBlocks>=0) = 1;
    
    % Hamming Decoder
    syndrome = mod((rBitBlock*H_transpose),2); % mod-2 airthmetic to get syndrome
    syndrome_Decimal = syndrome*([4 2 1]'); % Syndrome as decimal equivalent
    
    % Creating error_vector
    % if the sydrome decimal value matches the 1-7, then that position
    % would become 1, while the rest should stay false
    error_vector = error_index(syndrome_Decimal+1)==1:7;
    c_hat = mod(rBitBlock+error_vector,2); % estimated codeword
    u_hat = c_hat(:,4:7); % estimated message word
    n_errors(i) = sum(u(:)~=u_hat(:));
end
% Simulated BER
% N*k message bits
% N*n code bits
simulatedBER = n_errors/(N*k);
% Theoretical BER
BER_Theoretical = qfunc(sqrt(2*Eb_N0));
% Plotting
figure
semilogy(Eb_N0_dB,simulatedBER,'-');
hold on
semilogy(Eb_N0_dB,BER_Theoretical,'--');
legend('Simulation', 'Theoretical');
xlabel('Eb/N0 (dB)')
ylabel('Bit Error Rate (BER)')
title('BER vs Eb/N0 curve')
hold off
