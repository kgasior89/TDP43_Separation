%Code is used to generate phase separation simulation. File will run
%independently from user input or data files. Will output zipped .mat file
%for all species. Each page in the matrix for each species is a different
%point in time. X&Y on each page is the 2D snapshot at that moment. 
%
%Code was modified by Dr. Gasior and Dr. Newby for TDP-43 and HSP70 
%form anisotropic intranuclear liquid spherical annuli paper. 
%This code will generate simul1ations for with TDP43 phase
%separation. 
%Email: kgasior@fsu.edu

function  LLSolver()
% Control Parameters
Ny = 128;
Nz = 128;
StableC0 = 10.;

dt    = 1.0e-3;
t_end = 1000;
time  = 0.;
round = 0;

Ly = 1.;
Lz = 1.;

hy = Ly/Ny;
hz = Lz/Nz;

%Model parameters
eps2    = 1.0e-4;
lam_P2  = 1.25e-1;
lam_R   = 3e-2;
lam_P   = 2.5e-1;
lam_Y   = 1.25e-1;
c1 = 1e0;
c2 = 1e-2;

sig = 0.025;
A = 0.2;

chi_P2Y = 2.5;
chi_P2R = 4.25;
chi_YR = 4.25;
chi_P2S = 4.25;
chi_YS = 4.25;
chi_RS = 1;



%Set initial conditions
P2 = zeros(Ny,Nz);
P  = ((0.5).*rand(Ny,Nz));
Y  = ((0.5).*rand(Ny,Nz));
R = 1-(P+Y);


% set data for old time steps
P2_old = P2;
R_old = R;
P_old = P;
Y_old = Y;





%Use latter for FFT transform
N = Ny; M = Nz;
Leig  = (((2*cos(pi*(0:N-1)'/(N)))-2)*ones(1,M)) +(ones(N,1)*((2*cos(pi*(0:M-1)/(M)))-2));
Leigh= Leig/(hy*hz);
Seig = Leig/(hy*hz);


% 
% start loops
while (time < t_end)
    
    %update time and rounds
    time = time + dt;
    round= round +1;
    
    %extropolation data
    P2_bar = 2*P2 - P2_old;
    R_bar  = 2*R - R_old;
    P_bar  = 2*P - P_old;
    Y_bar  = 2*Y - Y_old;

%     
    %Solver go-go-go
    P2_np1 = solve_P2_fft(P2,P2_old,P2_bar,P_bar,R_bar,Y_bar);
    P_np1  = solve_P_fft(P,P_old,P_bar,P2_bar);
    R_np1 = solve_R_fft(R,R_old,R_bar,P2_bar,Y_bar);
    Y_np1 = solve_Y_fft(Y,Y_old,Y_bar,P2_bar,R_bar);
    
    %update data
    P2_old = P2; P2 = P2_np1;
    P_old = P;   P = P_np1;
    R_old = R;   R = R_np1;
    Y_old = Y;   Y = Y_np1;

%     
    if mod(round,100) == 0
        sprintf('Visualize time t=%f',round*dt)    
%         subplot(2,2,1);image(P,'CDataMapping','scaled');colormap(jet);colorbar;caxis([0 1]);xlabel('P');set(gca,'Xticklabel',[]);set(gca,'Yticklabel',[]);
%         subplot(2,2,2);image(R,'CDataMapping','scaled');colormap(jet);colorbar;caxis([0 1]);xlabel('R');set(gca,'Xticklabel',[]);set(gca,'Yticklabel',[]);
%         subplot(2,2,3);image(P2,'CDataMapping','scaled');colormap(jet);colorbar;caxis([0 1]);xlabel('P2');set(gca,'Xticklabel',[]);set(gca,'Yticklabel',[]);
%         subplot(2,2,4);image(Y,'CDataMapping','scaled');colormap(jet);colorbar;caxis([0 1]);xlabel('Y');set(gca,'Xticklabel',[]);set(gca,'Yticklabel',[]);
        page = round/100;
        P_mat(:,:,page) = P;
        P2_mat(:,:,page) = P2;
        R_mat(:,:,page) = R;
        Y_mat(:,:,page) = Y;
        pause(2);
    end

end

save('P.mat','P_mat');
zip('P.zip','P.mat');
delete P.mat
save('P2.mat','P2_mat');
zip('P2.zip','P2.mat');
delete P2.mat
save('R.mat','R_mat');
zip('R.zip','R.mat');
delete R.mat
save('Y.mat','Y_mat');
zip('Y.zip','Y.mat');
delete Y.mat

    function P2_np1 = solve_P2_fft(P2,P2_old,P2_bar,P_bar,R_bar,Y_bar)
        A_P2 = 1.5/dt * ones(N,M) - (lam_P2*StableC0*Leigh) + (lam_P2*eps2*Leigh.*Leigh);
        F_P2 = lam_P2* (chemical_potential_P2(P2_bar,Y_bar,R_bar)-StableC0*P2_bar);
        hat_rhs = 1./(2*dt)*(4*dct2(P2)-dct2(P2_old))+ Seig.*dct2(F_P2) + dct2(reactive_P2(P2_bar,P_bar));
        hat_P2 = hat_rhs./A_P2;
        P2_np1 = idct2(hat_P2);
    end

    function R_np1 = solve_R_fft(R,R_old,R_bar,P2_bar,Y_bar)
        A_R = 1.5/dt * ones(N,M) - (lam_R*StableC0*Leigh) + (lam_R*eps2*Leigh.*Leigh);
        F_R = lam_R* (chemical_potential_R(R_bar,P2_bar,Y_bar)-StableC0*R_bar);
        hat_rhs = 1./(2*dt)*(4*dct2(R)-dct2(R_old))+ Seig.*dct2(F_R);
        hat_R = hat_rhs./A_R;
        R_np1 = idct2(hat_R);
    end

    function Y_np1 = solve_Y_fft(Y,Y_old,Y_bar,P2_bar,R_bar)
        A_Y = 1.5/dt * ones(N,M) - (lam_Y*StableC0*Leigh) + (lam_Y*eps2*Leigh.*Leigh);
        F_Y = lam_Y* (chemical_potential_Y(Y_bar,P2_bar,R_bar)-StableC0*Y_bar);
        hat_rhs = 1./(2*dt)*(4*dct2(Y)-dct2(Y_old))+ Seig.*dct2(F_Y);
        hat_Y = hat_rhs./A_Y;
        Y_np1 = idct2(hat_Y);
    end

    function P_np1 = solve_P_fft(P,P_old,P_bar,P2_bar)
        %solver for Whi3
        A_P = 1.5/dt * ones(N,M) - lam_P * Leigh;
        hat_rhs = 1./(2*dt)*(4*dct2(P)-dct2(P_old)) + dct2(reactive_P(P2_bar,P_bar));
        hat_P = hat_rhs./A_P;
        P_np1 = idct2(hat_P);
    end

function Fr_P2 = chemical_potential_P2(x,Y_bar,R_bar)
    
       L1 = 1.*(x < sig);
       L2 = 1.*((1-(x+Y_bar+R_bar)) < sig);
       M1 = L1.*sig + (1-L1).*x;
       M2 = L2.*sig + (1-L2).*(1-(x+Y_bar+R_bar));
       R11 = L1.*(log(sig) + x./sig - log(M2) - 1); 
       R21 = L2.*(log(M1) + 1 - log(sig) - ((1-(x+Y_bar+R_bar))./sig));
       R31 = (1- L1).*(1-L2).*(log(M1) - log(M2));
       
       Fr_P2 = A.*(R11 + R21 + R31 + chi_P2Y*Y_bar + chi_P2R*R_bar - chi_P2S*x ...
           -chi_RS*R_bar-chi_YS*Y_bar+ chi_P2S*(1 - (x+Y_bar+R_bar)));
end

function Fr_R = chemical_potential_R(x,P2_bar,Y_bar)
    
       L12 = 1.*(x < sig);
       L22 = 1.*((1-(x+P2_bar+Y_bar)) < sig);
       M12 = L12.*sig + (1-L12).*x;
       M22 = L22.*sig + (1-L22).*(1-(x+P2_bar+Y_bar));
       R12 = L12.*(log(sig) + x./sig - log(M22) - 1); 
       R22 = L22.*(log(M12) + 1 - log(sig) - ((1-(x+P2_bar+Y_bar))./sig));
       R32 = (1- L12).*(1-L22).*(log(M12) - log(M22));
       
       Fr_R = A.*(R12 + R22 + R32 + chi_P2R*P2_bar + chi_YR*Y_bar - chi_P2S*P2_bar...
           -chi_YS*Y_bar-chi_RS*x+ chi_RS*(1 -(x+P2_bar+Y_bar)));
end

function Fr_Y = chemical_potential_Y(x,P2_bar,R_bar)

       L13 = 1.*(x < sig);
       L23 = 1.*((1-(x+P2_bar+R_bar)) < sig);
       M13 = L13.*sig + (1-L13).*x;
       M23 = L23.*sig + (1-L23).*(1-(x+P2_bar+R_bar));
       R13 = L13.*(log(sig) + x./sig - log(M23) - 1); 
       R23 = L23.*(log(M13) + 1 - log(sig) - ((1-(x+P2_bar+R_bar))./sig));
       R33 = (1- L13).*(1-L23).*(log(M13) - log(M23));
       
       Fr_Y = A.*(R13 + R23 + R33 + chi_P2Y*P2_bar + chi_YR*R_bar -chi_P2S*P2_bar...
           -chi_RS*R_bar - chi_YS*x + chi_YS*(1 - (x + P2_bar+R_bar)));
end


    function rhs_P2 = reactive_P2(P2_bar,P_bar)
        %reactive terms for N1
        rhs_P2 = c1*P_bar.*P_bar-c2*P2_bar;
    end


    function rhs_P = reactive_P(P2_bar,P_bar)
        %reactive terms for Whi3
        rhs_P = -2*c1*P_bar.*P_bar+2*c2*P2_bar;
    end


end


