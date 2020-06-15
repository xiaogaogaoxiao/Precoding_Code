clear,clc

%Nt = 144;
%Nr = 16;
%Ns = 4;
%NRF = 4;
Nt = 144;
Nr = 16;
Ns = 2;
NRF = 8;
V_RF=eye(Nt,NRF);
W_RF=eye(Nr,NRF);
SNR_dB = -20:2:0;
SNR = 10.^(SNR_dB./10);
realization = 2;
smax = length(SNR); % enable the parallel

Nc = 30;
c = 1/sqrt(Nc)*exp(1i*[0:2*pi/Nc:2*pi-2*pi/Nc])';
C = kron(eye(NRF),c);
P=sqrt(Nr*Nt);
%parfor reali = 1:realization
for reali = 1:realization
    reali
    [ H,At,Ar,Fopt,Wopt ] = channel_realization(Nt,Nr,Ns);
    
    [ FRFA, FBBA ] = AE_AltMin( Fopt, NRF);
    FBBA = sqrt(Ns) * FBBA / norm(FRFA * FBBA,'fro');
    [ WRFA, WBBA ] = AE_AltMin( Wopt, NRF);
    
    [ FRFO, FBBO ] = OMP( Fopt, NRF, At);
    FBBO = sqrt(Ns) * FBBO / norm(FRFO * FBBO,'fro');
    [ WRFO, WBBO ] = OMP( Wopt, NRF, Ar);

    [ FRF, FBB ] = my_AltMin_new( Fopt, C);
    FBB = sqrt(Ns) * FBB / norm(FRF * FBB,'fro');
    [ WRF, WBB ] = my_AltMin_new( Wopt, C);
    
    for s = 1:smax
%     [ V_t,FRF1, FBB1 ] = hybrid_beamforming_vt_w( H,NRF,V_RF,P,Nt,P/SNR(s));
%     
%     [ W_t] = hybrid_beamforming_wt_w(H,NRF,W_RF,V_t,Nr,P,P/SNR(s));
%     R1(s,reali)=log2(det(eye(Nr)+SNR(s)*W_t*((W_t'*W_t)^(-1))*W_t'*H*V_t*V_t'*H'));
    
        RA(s,reali) = log2(det(eye(Ns) + SNR(s) * pinv(WRFA*WBBA) * H * FRFA * FBBA * FBBA' * FRFA' * H' * WRFA*WBBA));
        ROM(s,reali) = log2(det(eye(Ns) + SNR(s) * pinv(WRFO*WBBO) * H * FRFO * FBBO * FBBO' * FRFO' * H' * WRFO*WBBO));
        R(s,reali) = log2(det(eye(Ns) + SNR(s) * pinv(WRF*WBB) * H * FRF * FBB * FBB' * FRF' * H' * WRF*WBB));
        RO(s,reali) = log2(det(eye(Ns) + SNR(s) * pinv(Wopt) * H * Fopt * Fopt' * H' * Wopt));
    end
end
plot(SNR_dB,sum(RO,2)/realization,'r-o','LineWidth',1.5);
hold on
grid on
plot(SNR_dB,sum(R,2)/realization,'m-^','LineWidth',1.5);
plot(SNR_dB,sum(RA,2)/realization,'Marker','>','Color',[0 0.498039215803146 0],'LineWidth',1.5);
plot(SNR_dB,sum(ROM,2)/realization,'b-v','LineWidth',1.5);
% plot(SNR_dB,sum(R1,2)/realization,'r-','LineWidth',1.5);
legend('Fully digital','Proposed FPS-AltMin','MO-AltMin [3]','OMP [2]')