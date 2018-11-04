N = 1e6;
sigma1_init = 2*rand(3,N)-ones(3,N);
sigma1_final = zeros(3,N);
sigma2_init =  2*rand(3,N)-ones(3,N);
sigma2_final = zeros(3,N);
axis_real = 2*rand(3,N)-ones(3,N);
ang_disp_real = pi*rand(1,N);
axis_arb = zeros(3,N);
ang_disp_arb = zeros(1,N);
ang2m_arb = zeros(1,N);
tic
for i = 1:N
    sigma1_init(:,i) = sigma1_init(:,i)./norm(sigma1_init(:,i));
    axis_real(:,i) = axis_real(:,i)./norm(axis_real(:,i));
    [sigma1_final(:,i),~] = Rotation3D(sigma1_init(:,i),axis_real(:,i),ang_disp_real(i));
    [sigma2_final(:,i),~] = Rotation3D(sigma2_init(:,i),axis_real(:,i),ang_disp_real(i));
    [axis_arb(:,i),ang_disp_arb(i),ang2m_arb(i)] = FindRotationAxis(sigma1_init(:,i),sigma1_final(:,i),sigma2_init(:,i),sigma2_final(:,i));
end
toc
d_axis = sqrt(sum((axis_arb-axis_real).^2,1));
d_ang = ang_disp_arb - ang_disp_real;
