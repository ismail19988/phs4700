%
% Devoir 3
% Billard fran�ais
% 
%
xtable=[0 2.8 2.8   0   0];
ytable=[0   0  1.53  1.53   0];
for Simulation=1:4
  if Simulation == 1
    disp('Simulation 1');
    xyb=[0.5;0.3];
    xyr=[1.5;1.1];
    xyj=[1.45;1.31];
    Vb0=[1.7;0.5];
  elseif Simulation == 2
    disp('Simulation 2');
    xyb=[0.5;0.3];
    xyr=[1.5;1.1];
    xyj=[1.45;1.31];
    Vb0=[1.5074;1.3145];
  elseif Simulation == 3
    disp('Simulation 3');
    xyb=[0.5;0.3];
    xyr=[1.5;1.1];
    xyj=[1.45;1.31];
    Vb0=[0.4;0.3];
  elseif Simulation == 4
    disp('Simulation 4');
    xyb=[0.5;0.3];
    xyr=[1.5;1.1];
    xyj=[1.45;1.31];
    Vb0=[1.6;1.35];
  end
  clf;
  hold;
  xlabel('x(m)')
  ylabel('y(m)')
  axis equal
  fill(xtable,ytable,'w');
 [coll tr tt xb yb xr yr xj yj] = Devoir3(xyb, xyr, xyj, Vb0);
 fprintf('Instant de transaction roulement avec-sans glissement \n');
  % disp(tr) 
  sz = size(tt,2);
  disp(xr(1:sz));
  plot(xb(1:sz),yb(1:sz),'bo',xr(1:sz),yr(1:sz),'bo',xj(1:sz),yj(1:sz),'bo')
  sz = size(tt,1);
%  plot(xyb(1),xyb(2),'kx')
%  fill(xyr(1),xyr(2),'kx')
%  fill(xyj(1),xyj(2),'kx')
  plot(xb(1:sz),yb(1:sz),'b')
  plot(xr(1:sz),yr(1:sz),'b')
  plot(xj(1:sz),yj(1:sz),'b')
  if coll == 0
    fprintf('Aucune bille touch�e \n');
  elseif coll == 1
    fprintf('Une seule bille touch�e \n');
  elseif coll == 2
    fprintf('Les deux billes sont touch�es \n');
    fprintf('Le point est marqu�\n');
  end
  hold;
  pause;
  clearvars xyb xyr xyj Vb0 coll tr tt xb yb xr yr xj yj
end;
