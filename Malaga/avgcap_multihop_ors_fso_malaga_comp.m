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
g = 2.7; % Pointing error coefficient 
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
    omega2 = 0.5*Cn2*(k_f^(7/6))*(L^(11/6)); % Rytov variance
    alpha = (exp(0.49*(omega2)*((1+0.56*(omega2^1.2))^(-7/6)))-1)^(-1); % Large scale scattering parameter
    beta = round((exp(0.51*(omega2)*((1+0.69*(omega2^1.2))^(-5/6)))-1)^(-1)); % Small scale scattering parameter
    b_o = 0.1079;
    rho = 0.596;
    omega = 1.3265;
    phi_A_phi_B = pi/2;
    y = 2*b_o*(1-rho);
    omega_1 = omega+(2*b_o*rho)+(2*cos(phi_A_phi_B)*sqrt(2*b_o*rho*omega));
    A_m = ((2*(alpha^(alpha/2)))*(((y*beta)/((y*beta)+omega_1))^(beta+(alpha/2))))/(gamma(alpha)*(y^(1+(alpha/2))));
    v = 4343/(L*psi);
    vv = (v^k)/((1+v)^k);
    chi = (A_m*alpha*beta*(g^2))/(2*A_o*((y*beta)+omega_1));
    zeta = (alpha*beta)/(A_o*((y*beta)+omega_1));
    B = A_o*gg*vv*(y+omega_1);
    for i = 1:length(Eb_N0_dB)
        cap_t(i) = 0;
        for d = 1:beta
            a_d = nchoosek(beta-1,(d)-1)*(((y*beta)+omega_1)^(1-((d)/2)))*(((omega_1)/y)^((d)-1))*((alpha/beta)^((d)/2))/factorial((d)-1);
            b_d = a_d*((((y*beta)+omega_1)/(alpha*beta))^((alpha+d)/2));
            const_1 = (2^(alpha+d-k-3))/(pi*log(2));
            const_2 = ((v^(k))*chi)/zeta;
            arg(i) = ((zeta^2)*(B^2))/(16*s_r*Eb_N0(i));
            cap_t(i) = cap_t(i)+(const_1*const_2*b_d*meijerG([0],[1,(1+g^2)/2,(2+g^2)/2,(1+v)/2,(1+v)/2,(2+v)/2,(2+v)/2],[(g^2)/2,(1+g^2)/2,alpha/2,(1+alpha)/2,d/2,(1+d)/2,v/2,v/2,(1+v)/2,(1+v)/2,0,0],[],arg(i)));
        end
    end
    cap_t_final(m,:) = cap_t;
end

%==============================Simulation=================================%

Nhop = [2,4,6];
for m = 1:3
    N = Nhop(m);
    L1 = 2000/(2*N); % Link distance form transmitter to IRS
    L2 = L1; % Link distance form IRS to reciever
    L = L1+L2; % Total link distance
    Cn2 = 2.3*10^(-13); % Refractive index parameter
    lambda = 1.55*10^(-6); % Wavelength
    k_f = (2*pi)/lambda; % Optical wavenumber
    omega2 = 0.5*Cn2*(k_f^(7/6))*(L^(11/6)); % Rytov variance
    alpha = (exp(0.49*(omega2)*((1+0.56*(omega2^1.2))^(-7/6)))-1)^(-1); % Large scale scattering parameter
    beta = round((exp(0.51*(omega2)*((1+0.69*(omega2^1.2))^(-5/6)))-1)^(-1)); % Small scale scattering parameter
    b_o = 0.1079;
    rho = 0.596;
    omega = 1.3265;
    phi_A_phi_B = pi/2;
    y = 2*b_o*(1-rho);
    omega_1 = omega+(2*b_o*rho)+(2*cos(phi_A_phi_B)*sqrt(2*b_o*rho*omega));
    A_m = ((2*(alpha^(alpha/2)))*(((y*beta)/((y*beta)+omega_1))^(beta+(alpha/2))))/(gamma(alpha)*(y^(1+(alpha/2))));
    v = 4343/(L*psi);
    vv = (v^k)/((1+v)^k);
    B = A_o*gg*vv*(y+omega_1);
    for i = 1:length(Eb_N0_dB)
        for j = 1:N
            sigma = sqrt(b_o*(1-rho));
            A1 = sigma.*(randn(1,M)+1i*randn(1,M));
            Z = sqrt(gamrnd(beta,omega_1/beta,1,M));
            G = gamrnd(alpha,1/alpha,1,M);
            h_a = ((abs(A1+Z)).^2).*G; % Atmospheric turbulance (Malaga distribution)
            r_d = raylrnd(sigma_j,1,M); % Radial displacement (Rayleigh distribution) 
            h_pe = A_o.*exp(-(2.*(r_d.^2))./(w_zeq^2)); % Pointing error
            tau = gamrnd(k,psi,1,M); % Attenuation coefficient (Gamma distribution)
            h_f = exp(-(tau.*L)./4343); % Fog
            h = h_a.*h_pe.*h_f; % Overall channel coefficient
            inst_snr(j,:) = (Eb_N0(i).*(h.^2))./(B^2); % Instantaneous SNR
            cap_s(j,:) = (log(1+(s_r.*inst_snr(j,:)))./log(2));
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
semilogy(Eb_N0_dB,cap_s_oall_final(1,:),'x',Eb_N0_dB,cap_s_oall_final(2,:),'v',Eb_N0_dB,cap_s_oall_final(3,:),'pentagram','linewidth',3,'markersize',8);
hold on
semilogy(Eb_N0_dB,cap_t_final,'k-','lineWidth',1.5);
legend('N = 2','N = 4','N = 6','Theoretical','Location','southeast')
grid on
xlabel('Average SNR (dB)')
ylabel('Ergodic Capacity')