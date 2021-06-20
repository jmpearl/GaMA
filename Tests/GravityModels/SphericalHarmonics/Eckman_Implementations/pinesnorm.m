

function accel = pinesnorm(MU, REQ, x, CNM, SNM, NMAX, MMAX, rnp)
R_F = rnp*x; %RAE
RMAG = norm(R_F);
S = R_F(1)/RMAG;
T = R_F(2)/RMAG;
U = R_F(3)/RMAG;
ANM = zeros(NMAX+3,NMAX+3); %RAE
ANM(1,1) = sqrt(2.0D0); %norm
for M = 0:NMAX+2 %RAE
    M_A = M + 1;
    if(M ~= 0) % DIAGONAL RECURSION
        ANM(M_A,M_A) = sqrt(1+(1/(2*M)))*ANM(M_A-1,M_A-1); %norm
    end
    if(M ~= NMAX+2) % FIRST OFF-DIAGONAL RECURSION %RAE
        ANM(M_A+1,M_A) = sqrt(2*M+3)*U*ANM(M_A,M_A); %norm
    end
    if(M < NMAX+1) % COLUMN RECURSION %RAE
        for N = M+2:NMAX+2 %RAE
            N_A = N + 1;
            ALPHA_NUM = (2*N+1)*(2*N-1);
            ALPHA_DEN = (N-M)*(N+M);
            ALPHA = sqrt(ALPHA_NUM/ALPHA_DEN);
            BETA_NUM = (2*N+1)*(N-M-1)*(N+M-1);
            BETA_DEN = (2*N-3)*(N+M)*(N-M);
            BETA = sqrt(BETA_NUM/BETA_DEN);
            ANM(N_A,M_A) = ALPHA*U*ANM(N_A-1,M_A) - BETA*ANM(N_A-2,M_A); %norm 
        end
    end
end
for N = 0:NMAX+2 %RAE
    N_A = N + 1;
    ANM(N_A,1) = ANM(N_A,1)*sqrt(0.5D0); %norm
end
RM = zeros(MMAX+2,1); %RAE
IM = zeros(MMAX+2,1); %RAE
RM(1) = 0.0D0;
IM(1) = 0.0D0;
RM(2) = 1.0D0;
IM(2) = 0.0D0;
for M = 1:MMAX
    M_RI = M+2;
    RM(M_RI) = S*RM(M_RI-1) - T*IM(M_RI-1);
    IM(M_RI) = S*IM(M_RI-1) + T*RM(M_RI-1);
end

RHO = (MU)/(REQ*RMAG);
RHOP = (REQ)/(RMAG);
G1 = 0.0D0;
G2 = 0.0D0;
G3 = 0.0D0;
G4 = 0.0D0;
for N = 0:NMAX
    N_A = N + 1;
    G1TEMP =0.0D0;
    G2TEMP =0.0D0;
    G3TEMP =0.0D0;
    G4TEMP =0.0D0;
    SM = 0.5D0;
    
    if (N>MMAX) %RAE
        nmodel=MMAX; %RAE
    else %RAE
        nmodel=N; %RAE
    end %RAE
    for M = 0:nmodel %RAE
        M_A= M + 1;
        M_RI= M + 2; %RAE
        DNM = CNM(N_A,M_A)*RM(M_RI)+ SNM(N_A,M_A)*IM(M_RI);
        ENM = CNM(N_A,M_A)*RM(M_RI-1) + SNM(N_A,M_A)*IM(M_RI-1);
        FNM = SNM(N_A,M_A)*RM(M_RI-1) - CNM(N_A,M_A)*IM(M_RI-1);
        ALPHA = sqrt(SM*(N-M)*(N+M+1)); %norm
        G1TEMP = G1TEMP + ANM(N_A,M_A)*(M)*ENM;
        G2TEMP = G2TEMP + ANM(N_A,M_A)*(M)*FNM;
        G3TEMP = G3TEMP + ALPHA*ANM(N_A,M_A+1)*DNM; %norm
        G4TEMP = G4TEMP + ((N+M+1)*ANM(N_A,M_A)+ALPHA*U*ANM(N_A,M_A+1))*DNM;%norm
        if(M == 0)
            SM = 1.0D0;
        end %norm
    end
    RHO=RHOP*RHO;
    G1=G1 + RHO*G1TEMP;
    G2=G2 + RHO*G2TEMP;
    G3=G3 + RHO*G3TEMP;
    G4=G4 + RHO*G4TEMP;
 
end
G_F = [G1 - G4*S;G2 - G4*T;G3 - G4*U];
accel=rnp'*G_F; %RAE
return
