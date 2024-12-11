clc;
clear all;

Eb_N0_dB = 0:5:40; % Avg SNR dB
Eb_N0 = 10.^(Eb_N0_dB./10); % Avg SNR
M = 10^6; % No. of bits or symbols
a_z = 0.1; % Aperture radius
w_z = 2.5; % Beamwidth
d = (sqrt(pi)*a_z)/(sqrt(2)*w_z);
w_zeq = sqrt(((w_z^2)*sqrt(pi)*erf(d))/(2*d*exp(-(d^2)))); % Equivalent beamwidth
A_o = (abs(erf(d)))^2; % Amount of power collected at r_d = 0
g = 2.7; %  pointing error coefficient 
gg = (g^2)/((g^2)+1);
sigma_j = w_zeq/(2*g); % Jitter standard deviation
k = 2; % Shape parameter of fog
psi = 13.12; % Scale parameter of fog
s_r = exp(1)/(2*pi);

%==============================Theoretical================================%

Nhop = [2,4,6];
for m = 1:3
    N = Nhop(m);
    L1 = 2000/(2*N); % Link distance form transmitter to IRS
    L2 = L1; % Link distance form IRS to reciever
    L = L1+L2; % Total link distance
    Cn2 = 2.3*10^(-13); % Refractive index parameter
    lambda = 1.55*10^(-6); % Wavelength
    k_f = (2*pi)/lambda; % Optical wavenumber
    omega2 = 0.5*Cn2*(k_f^(7/6))*(L1^(11/6)); % Rytov variance
    alpha = (exp(0.49*(omega2)*((1+0.56*(omega2^1.2))^(-7/6)))-1)^(-1); % Large scale scattering parameter
    beta = (exp(0.51*(omega2)*((1+0.69*(omega2^1.2))^(-5/6)))-1)^(-1); % Small scale scattering parameter
    lambda_a = alpha+2;
    A = 1/(gamma(alpha)*gamma(beta)*gamma(lambda_a));
    v = 4343/(L1*psi);
    vv = (v^k)/((1+v)^k);
    chi = (A*alpha*beta*(g^2))/(A_o*(1-lambda_a));
    zeta = (alpha*beta)/(A_o*(1-lambda_a));
    B = (A_o*gg*vv)^2; % Expectation value of overall channel coefficient
    for i = 1:length(Eb_N0_dB)
        const_1 = (2^((2*alpha)+(2*beta)+(2*lambda_a)-(2*k)-8))/((pi^3)*log(2));
        const_2 = ((v^(2*k))*(chi^2))/(zeta^2);
        arg(i) = ((zeta^4)*(B^2))/(16*s_r*Eb_N0(i));
        cap_t(i) = const_1*const_2*meijerG([(1-lambda_a)/2,(2-lambda_a)/2,(1-lambda_a)/2,(2-lambda_a)/2,0],[1,(1+g^2)/2,(2+g^2)/2,(1+g^2)/2,(2+g^2)/2,(1+v)/2,(1+v)/2,(2+v)/2,(2+v)/2,(1+v)/2,(1+v)/2,(2+v)/2,(2+v)/2],[(g^2)/2,(1+g^2)/2,(g^2)/2,(1+g^2)/2,alpha/2,(1+alpha)/2,alpha/2,(1+alpha)/2,beta/2,(1+beta)/2,beta/2,(1+beta)/2,v/2,v/2,(1+v)/2,(1+v)/2,v/2,v/2,(1+v)/2,(1+v)/2,0,0],[],arg(i));
    end
    cap_t_final(m,:) = cap_t;
end

%==============================Simulation=================================%

for m = 1:3
    N = Nhop(m);
    L1 = 2000/(2*N); % Link distance form transmitter to IRS
    L2 = L1; % Link distance form IRS to reciever
    L = L1+L2; % Total link distance
    Cn2 = 2.3*10^(-13); % Refractive index parameter
    lambda = 1.55*10^(-6); % Wavelength
    k_f = (2*pi)/lambda; % Optical wavenumber
    omega2 = 0.5*Cn2*(k_f^(7/6))*(L1^(11/6)); % Rytov variance
    alpha = (exp(0.49*(omega2)*((1+0.56*(omega2^1.2))^(-7/6)))-1)^(-1); % Large scale scattering parameter
    beta = (exp(0.51*(omega2)*((1+0.69*(omega2^1.2))^(-5/6)))-1)^(-1); % Small scale scattering parameter
    lambda_a = alpha+2;
    A = 1/(gamma(alpha)*gamma(beta)*gamma(lambda_a));
    v = 4343/(L1*psi);
    vv = (v^k)/((1+v)^k);
    B = (A_o*gg*vv)^2; % Expectation value of overall channel coefficient
    for i = 1:length(Eb_N0_dB)
        for j = 1:N
            h_a = gamrnd(alpha,1/alpha,2,M).*gamrnd(beta,1/beta,2,M).*(1./gamrnd(lambda_a,1/(lambda_a-1),2,M)); % Atmospheric turbulance (IGGG distribution)
            r_d = raylrnd(sigma_j,2,M); % Radial displacement (Rayleigh distribution) 
            h_pe = A_o.*exp(-(2.*(r_d.^2))./(w_zeq^2)); % Pointing error
            tau = gamrnd(k,psi,2,M); % Attenuation coefficient (Gamma distribution)
            h_f = exp(-(tau.*L1)./4343); % Fog
            h_i = h_a.*h_pe.*h_f; % Channel coefficient from transmitter to IRS
            h(j,:) = prod(h_i); % Overall channel coefficient
            inst_snr(j,:) = (Eb_N0(i).*(h(j,:).^2))./(B^2); % Instantaneous SNR
            cap_s(j,:) = (log(1+(s_r.*inst_snr(j,:)))./(log(2)));
            cap_s_mc(j) = sum(cap_s(j,:))./M;
        end
        if N == 1
            cap_s_oall(i,:) = cap_s_mc;
        else
            cap_s_oall(i,:) = min(cap_s_mc);
        end
    end
    cap_s_oall_final(m,:) = cap_s_oall;
end

figure
semilogy(Eb_N0_dB,cap_s_oall_final(1,:),'^',Eb_N0_dB,cap_s_oall_final(2,:),'square',Eb_N0_dB,cap_s_oall_final(3,:),'o','linewidth',3,'markersize',8);
hold on
semilogy(Eb_N0_dB,cap_t_final,'k-','lineWidth',1.5);
legend('N = 2','N = 4','N = 6','Theoretical','Location','southeast')
grid on
xlabel('Average SNR (dB)')
ylabel('Ergodic Capacity')