% Examine the  eigenvalues of the cable sail matrix and plot
% Written by:   Keegan Bunker 
%               PhD student, University of Minnesota 
%               bunke029@umn.edu
%
% Last major modifications: May 31, 2023

%LoadConstants

A1 = zeros(size(const.M));
A2 = diag(ones(1,length(const.M)));
A3 = -const.invM * const.K;
A4 = -const.invM*const.D;

Adamp = [A1 A2; A3 A4];

A = [A1 A2; A3 A1];

eigsdmp = eig(Adamp)';
eigs = eig(A)';

figure()
plot(eigsdmp, 'ro')
hold on
plot(eigs, 'bo')
grid on 
legend('Damped', 'Undamped','Location', 'NW','Interpreter', 'latex','FontSize',16)


title('Eigenvalues','Interpreter', 'latex','FontSize',18)
xlabel('Real','Interpreter', 'latex','FontSize',16)
ylabel('Imaginary','Interpreter', 'latex','FontSize',16)

% list the angles of eigenvalues
angledmp = angle(eigsdmp);
anglesundmp = angle(eigs);
allangles = [angledmp; anglesundmp; angledmp - anglesundmp];



disp('eigen values, damped and undamped')
[eigsdmp',eigs']

disp('Damped eigenvalue angles')
disp('Undamped eigenvalue angles')
disp('Difference in angles')
allangles
