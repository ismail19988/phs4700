% Ismail Bakkouri - 1954157
% Mohammed Ariful Islam - 1950221
% John Maliha - 1984959
% Dawut Esse - 1956802

% NOTE : POUR EXECUTER LE DEVOIR, IL FAUT CHANGER LA LIGNE 46 DU ROULEDEVOIR2.M 
% REMPLACEZ
% sz=size(tt,2);
% PAR
% sz=size(tt,1);

function [coup Vbf t x y z] = Devoir2(xy0, Vb0, Wb0)

 % Constantes
Mb = 0.450;
Rb = 0.11;
A  = pi * (Rb ^ 2);
p = 1.2754;
mu = 1.8 * (10 ^ -5);
dt = 0.001;
t0 = 0;
coup = -1;
longueurTerrain = 120;
largeurTerrain = 90;
largeurBut = 7.32;
hauteutBut = 2.44;
erreur = 0.001;

% petit y
poteau1 = (largeurTerrain - largeurBut)/2;

% grand y
poteau2 = (largeurTerrain + largeurBut)/2;

% poteau horizontal
poteau3 = hauteutBut;

% calculer la force resultate appliquée sur la balle

    function [Force] = CalculerForce(vitesse)
        
        % force gravitationelle
        Force = transpose([0, 0, -Mb * 9.8]);
        
        % Force de frottement visqueux
        Re = p * norm(vitesse) * Rb / mu;
        
        if(Re < 100000)
            Cvis = 0.235 * norm(vitesse);
        elseif(100000 < Re && Re < 135000)
            Cvis = (0.235 * norm(vitesse)) - (0.125 * norm(vitesse) * ((Re - 100000)/35000));        
        else
            Cvis = 0.110 * norm(vitesse);        
        end
        
        Force = Force + ((-A) * p * Cvis * vitesse);
        
        % Force de Magnus (vitesse angulaire constante, donc depend seulement de la vitesse de translation)
        coefficient_magnus = 0.1925 * (((norm(Wb0) * Rb) / norm(vitesse)) ^ (1/4));
        Force = Force + (p * coefficient_magnus * A * (norm(vitesse)^2) * ((cross(Wb0, vitesse))/norm((cross(Wb0, vitesse))))); 
       
    end

    function [g] = g(q0, t0)
        acceleration = CalculerForce([q0(1); q0(2); q0(3);]) / Mb;
        g = [acceleration(1); acceleration(2); acceleration(3); q0(1); q0(2); q0(3);];
    end

% Vecteur de position X
x = [xy0(1)];

% Vecteur de position y
y = [xy0(2)];

% Vecteur de position Z
z = [0];

% Vecteur de temps
t = [0];

    function [done] = fini(position)
        
        % but 
        if((position(1) < 0 || position(1) > longueurTerrain) && ...
           (position(2) > poteau1 && position(2) < poteau2 && ...
           (position(3) + Rb < poteau3 && position(3) - Rb > 0)))
               coup = 0;
               done = true;
        % sol
        elseif (position(1) > 0 && position(1) < longueurTerrain && ...
            position(2) > 0 && position(2) < largeurTerrain ...
            && position(3) - Rb < 0)
            coup = 1;
            done = true;
            
        % poteau    
        elseif(((position(1) - Rb <= erreur || (position(1) + Rb >= longueurTerrain + erreur))) ...
             && ( ... 
                    (position(3) - Rb < (poteau3 - erreur) && position(3) + Rb > poteau3 + erreur) || ...
                    (position(2) - Rb < (poteau2 - erreur) && position(2) + Rb > poteau2 + erreur) || ...
                    (position(2) - Rb < (poteau1 - erreur) && position(2) + Rb > poteau1 + erreur)) ...
                )
               coup = 3;
               done = true;
               
         % sort du terrain         
         elseif (position(1) < 0 || position(1) > longueurTerrain) || (position(2) < 0 || position(2) > largeurTerrain)             
            coup = 2;
            done = true;
       
        % sinon, il faut continuer le calcul
        else 
            done = false;
            
        end
        
    end

% Calcul methode de Runge-Kutta
q0 = [Vb0(1); Vb0(2); Vb0(3); xy0(1); xy0(2); Rb];
qs = q0;

while true
    q0 = qs;
    position = [q0(4); q0(5); q0(6)];
    
    % on ajoute la position calculée dans nos vecteurs de sortie
    x =  [x; position(1);];
    y =  [y; position(2);];
    z =  [z; position(3);];
    t =  [t; t0;];
    
    % On verifie si la position recalculée permet de conclure le coup
    if(fini(position))
        break;
    end
    
    % Réévaluation des paramètres.
    k1 = g(q0, t0);
    k2 = g((q0 + (k1 * dt/2)), (t0 + dt/2));
    k3 = g((q0 + (k2 * dt/2)), (t0 + dt/2));
    k4 = g((q0 + (k3 * dt)), (t0 + dt));
    qs = q0 + dt / 6 * (k1 + (2 * k2) + (2 * k3) + k4);
    t0 = t0 + dt;
    
end

% Vecteur de vitesse finale de la balle.
Vbf = [qs(1); qs(2); qs(3)];


end