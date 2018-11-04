function [axis_out,ang_disp_out,ang2m] = FindRotationAxis(sigma1_init,sigma1_final,sigma2_init,sigma2_final)

%axis_0 = [0 0 1];
%ang_disp_0 = 0.5;
%max_ang_disp = 0.1; % maximum angular displacement
%r0 = [axis_0 ang_disp_0];
%lb = [-1 -1 -1 0];
%ub = [1 1 1 max_ang_disp];
%numstart = 20;
%ms = MultiStart('StartPointsToRun','bounds-ineqs','UseParallel','always','Display','off');
%gs = GlobalSearch('StartPointsToRun','bounds-ineqs','Display','off');
%hybridopts = optimset('Algorithm','interior-point','TolFun',1e-4);
% confuneq can be adjusted to accounted for upper/lower bounds on
% displacement angle
%problem = createOptimProblem('fmincon','objective',@(r) find_ang2(sigma1_init,sigma1_final,sigma2_init,sigma2_final,r),'x0',r0,'lb',lb,'ub',ub,'nonlcon',@confuneq,'options',hybridopts);
%[r_out,ang2m] = run(ms,problem,numstart);
%[r_out,ang2m] = run(gs,problem);

P = [sigma1_init sigma2_init];
Q = [sigma1_final sigma2_final];

[R,ang2m] = Kabsch_mjh(P,Q);
e_vec = EulerParams(R);

ang_disp_out = 2*acos(e_vec(1));
axis_out = e_vec(2:4)/sin(ang_disp_out/2);