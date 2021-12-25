clear
%
% Procedure servant à rouler le devoir 1
%
% Cas 1
AngRotCas1=0.0;
vangulaireCas1=[0;0;0];
forceCas1=[5e6;0;0];
posCas1=[0;0;0];
[pcmCas1 ICas1 alphaCas1]=Devoir1(posCas1,AngRotCas1,vangulaireCas1,forceCas1);
fprintf('\nResultats missile verticale \n');
fprintf('  Centre de masse = ( %10.5f  %10.5f  %10.5f )\n',pcmCas1(1),pcmCas1(2),pcmCas1(3));
fprintf('  Moment inertie  =\n   %14.0f   %14.0f   %14.0f \\\\ \n   %14.0f   %14.0f   %14.0f \\\\ \n   %14.0f   %14.0f   %14.0f \\\\ \n',...
     ICas1(1,1),ICas1(1,2),ICas1(1,3),ICas1(2,1),ICas1(2,2),ICas1(2,3),ICas1(3,1),ICas1(3,2),ICas1(3,3));
fprintf('  acc angulaire   = ( %10.5f  %10.5f  %10.5f )\n',alphaCas1(1),alphaCas1(2),alphaCas1(3));
%
% Cas 2
AngRotCas2=0.5;
vangulaireCas2=[0.05;0;0.001];
forceCas2=[5e6;0.1;1.4];
posCas2=[0;-20.0;300];
[pcmCas2 ICas2 alphaCas2]=Devoir1(posCas2,AngRotCas2,vangulaireCas2,forceCas2);
fprintf('\nResultats missile incline\n');
fprintf('  Centre de masse = ( %10.5f  %10.5f  %10.5f )\n',pcmCas2(1),pcmCas2(2),pcmCas2(3));
fprintf('  Moment inertie  =\n   %14.0f   %14.0f  %14.0f \\\\ \n   %14.0f   %14.0f   %14.0f \\\\ \n   %14.0f   %14.0f   %14.0f \\\\ \n',...
     ICas2(1,1),ICas2(1,2),ICas2(1,3),ICas2(2,1),ICas2(2,2),ICas2(2,3),ICas2(3,1),ICas2(3,2),ICas2(3,3));
fprintf('  acc angulaire   = ( %10.5f  %10.5f  %10.5f )\n',alphaCas2(1),alphaCas2(2),alphaCas2(3));
