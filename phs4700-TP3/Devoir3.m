% Ismail Bakkouri - 1954157
% Mohammed Ariful Islam - 1950221
% John Maliha - 1984959
% Dawut Esse - 1956802

% NOTE : POUR EXECUTER LE DEVOIR, IL FAUT CHANGER LA LIGNE 46 DU ROULEDEVOIR3.M 
% REMPLACEZ
%  sz=size(tt,2);
%  % plot(xb(1:sz),yb(1:sz),'w',xr(1:sz),yr(1:sz),'r',xj(1:sz),yj(1:sz),'y')
% PAR
%  sz = size(tt,2);
%  disp(xr(1:sz));
%  plot(xb(1:sz),yb(1:sz),'wo',xr(1:sz),yr(1:sz),'ro',xj(1:sz),yj(1:sz),'yo')
%  sz = size(tt,1);


function [coll tr t xb yb xr yr xj yj]= Devoir3(xyb, xyr, xyj, Vb0)
    m = 0.210;
    r = 0.031;
    grav = 9.8;
    mu_r = 0.03;
    mu_c = 0.3;
    mu_s_mur = 0.4;
    mu_s_bille = 0.15;
    mu_c_bille = 0.1;
    dt = 0.0001;
    r_pc = [0; 0; -r];
    
    N = [0; 0; m * grav];
    I = [(2/5 * m * r^2) 0 0; 0 (2/5 * m * r^2) 0; 0 0 (2/5 * m * r^2);];
    restitution_bille_bord = 0.8;
    restitution_bille_bille = 0.9;
    
    hauteur_table = 1.53;
    largeur_table = 2.80;
    
    col_rouge = 0;
    col_jaune = 0;
    
    function [F_f, F_r] = calculerForces(vitesse_cm, glisse, v_roulement)
        if(glisse)
            v_glissement = vitesse_cm - v_roulement;
            F_f = -mu_c * norm(N) * v_glissement / norm(v_glissement);
            if(sum(isnan(F_f(:))))
                F_f = [0 0 0];
            end
            F_r = transpose([0 0 0]);
        else
            if(norm(vitesse_cm) > 0.001)
                F_r = -mu_r * norm(N) * vitesse_cm / norm(vitesse_cm);
            else
                F_r = transpose([0 0 0]);
            end
            F_f = transpose([0 0 0]);
        end
    end
    
    function [moment]= calculerMoments(vitesse, glisse, v_roulement, v_angulaire)
        [force_frottement, force_roulement] = calculerForces(vitesse, glisse, v_roulement);
        if(glisse)
            T_g = cross(r_pc, force_frottement);
            T_r = transpose([0 0 0]);
        else
            T_r = transpose([0 0 0]);
            T_g = transpose([0 0 0]);
        end
        
        moment = T_g + T_r;
    end

    function [g] = g(q0, glisse, v_roulement, v_angulaire)
        vitesse = [q0(1); q0(2); q0(3);];
        [force_frottement, force_roulement] = calculerForces(vitesse, glisse, v_roulement);
        forces = force_frottement + force_roulement;
        acc = forces / m;
        moment = calculerMoments(vitesse, glisse, v_roulement, v_angulaire);
        aa = inv(I) * moment;
        g = [acc(1); acc(2); acc(3); q0(1);  q0(2);  q0(3); aa(1); aa(2); aa(3);];
    end

    function [qs] = RK(q0, glisse, v_roulement, v_angulaire)
        k1 = g(q0, glisse, v_roulement, v_angulaire);
        k2 = g((q0 + (k1 * dt/2)), glisse, v_roulement, v_angulaire);
        k3 = g((q0 + (k2 * dt/2)), glisse, v_roulement, v_angulaire);
        k4 = g((q0 + (k3 * dt)), glisse, v_roulement, v_angulaire);
        qs = q0 + dt / 6 * (k1 + (2 * k2) + (2 * k3) + k4);
    end
    
    function [q0, blockers, glisse] = checkCollisions(q0_initial, blockers, glisse)
        q0 = q0_initial;
        vitesse_avant_collision = [q0_initial(1); q0_initial(2); q0_initial(3)];
        pos = [q0(4); q0(5); q0(6)];
        
        if(glisse)
            vitesse_angulaire_avant_collision = [q0(7); q0(8); q0(9)];
        else
            vitesse_angulaire_avant_collision = -cross(r_pc, vitesse_avant_collision)/(r^2);
        end
        
        collision = false; 
        if(pos(1) + r >= largeur_table && ~blockers(2))
            n = [-1; 0; 0];
            r_ap = [r; 0; 0];
            collision = true; 
            blockers(2) = true;
        elseif (pos(1) -  r <= 0 && ~blockers(4))
            n = [1; 0; 0];
            r_ap = [-r; 0; 0];
            collision = true; 
            blockers(4) = true;
       elseif (pos(2) +  r >= hauteur_table && ~blockers(1))
            n = [0; -1; 0];
            r_ap = [0; r; 0];
            collision = true;       
            blockers(1) = true;
        elseif (pos(2) -  r <= 0 && ~blockers(3))
            n = [0; -1; 0];
            r_ap = [0; -r; 0];
            collision = true; 
            blockers(3) = true;
        end
        
        if (collision)
            u = cross(vitesse_avant_collision, n)/norm(cross(vitesse_avant_collision, n));
            
            t_chapeau = cross(n, u);
            
            j = -m * (1 + restitution_bille_bord) * dot(n, vitesse_avant_collision);
            
            G_a_t = dot(t_chapeau, (inv(I) * cross(cross(r_ap, t_chapeau), r_ap)));
            alpha_t = 1/((1/m) + G_a_t);
            
            if (mu_s_mur * (1 + restitution_bille_bord) * abs(dot(n, vitesse_avant_collision)) < abs(dot(t_chapeau, vitesse_avant_collision)))
                jt = alpha_t * mu_c * (1 + restitution_bille_bord) * dot(n, vitesse_avant_collision);
            else
                jt = -alpha_t * abs(dot(t_chapeau, vitesse_avant_collision));
            end
            J = n * j + t_chapeau * jt;

            vitesse_apres_collision = vitesse_avant_collision + (J / m);

            vitesse_angulaire_apres_collision = vitesse_angulaire_avant_collision + (inv(I) * cross(r_ap, J));
            
            % mettre a jour 4 5 6
            q0 = updateQ0ApresCollision(q0, vitesse_apres_collision, vitesse_angulaire_apres_collision);
            vitesse_roulement = -cross(vitesse_angulaire_apres_collision, r_pc)
    
          if ((abs(norm(vitesse_apres_collision - vitesse_roulement))) > 0.0005)
              % glisse = true;
          end
            
        end
    end 
    
    function [J] = calculerImpulsion(vr_avant, n, r_ap, r_bp)
        u = cross(vr_avant, n)/norm(cross(vr_avant, n));
        t_chapeau = cross(n, u);
        
        G_a = dot(n, (inv(I) * cross(cross(r_ap, n), r_ap)));
        G_b = dot(n, (inv(I) * cross(cross(r_bp, n), r_bp)));
        alpha = 1/((1/m) + (1/m) + G_a + G_b);
        j = -alpha * (1 + restitution_bille_bille) * dot(n, vr_avant);
        
        G_a_t = dot(t_chapeau, (inv(I) * cross(cross(r_ap, t_chapeau), r_ap)));
        G_b_t = dot(t_chapeau, (inv(I) * cross(cross(r_bp, t_chapeau), r_bp)));
        alpha_t = 1/((1/m) + (1/m) + G_a_t + G_b_t);
          
        if (mu_s_bille * (1 + restitution_bille_bille) * abs(dot(n, vr_avant)) < abs(dot(t_chapeau, vr_avant)))
            jt = alpha_t * mu_c_bille * (1 + restitution_bille_bille) * dot(n, vr_avant);
        else
            jt = -alpha_t * abs(dot(t_chapeau, vr_avant));
        end
        J = n * j + t_chapeau * jt;
    end

    function [q0_blanche, q0_rouge, q0_jaune] = checkCollisionsBille(q0_blanche_initial, q0_rouge_initial, q0_jaune_initial)
        q0_blanche = q0_blanche_initial;
        q0_rouge = q0_rouge_initial;
        q0_jaune = q0_jaune_initial;
        
        position_blanche = [q0_blanche(4); q0_blanche(5); q0_blanche(6)];
        position_rouge   = [q0_rouge(4);   q0_rouge(5); q0_rouge(6)];
        position_jaune   = [q0_jaune(4);   q0_jaune(5); q0_jaune(6)];
        vitesse_avant_collision_blanche = [q0_blanche(1); q0_blanche(2); q0_blanche(3)];
        vitesse_avant_collision_rouge = [q0_rouge(1); q0_rouge(2); q0_rouge(3)];
        vitesse_avant_collision_jaune = [q0_jaune(1); q0_jaune(2); q0_jaune(3)];
        
        if(blanche_glisse)
            vitesse_angulaire_avant_collision_blanche = [q0_blanche(7); q0_blanche(8); q0_blanche(9)];
        else
            vitesse_angulaire_avant_collision_blanche = -cross(r_pc, vitesse_avant_collision_blanche)/(r^2);
        end
        
        if(rouge_glisse)
            vitesse_angulaire_avant_collision_rouge = [q0_rouge(7); q0_rouge(8); q0_rouge(9)];
        else
            vitesse_angulaire_avant_collision_rouge = -cross(r_pc, vitesse_avant_collision_rouge)/(r^2);
        end
        
        if(jaune_glisse)
            vitesse_angulaire_avant_collision_jaune = [q0_jaune(7); q0_jaune(8); q0_jaune(9)];
        else
            vitesse_angulaire_avant_collision_jaune = -cross(r_pc, vitesse_avant_collision_jaune)/(r^2);
        end
        
        if (norm(position_blanche - position_rouge) <= r + r && ~block_blanche_rouge)
          blockers_blanc = resetBlockers();
          blockers_rouge = resetBlockers();
          block_blanche_rouge = true;
          n = (position_blanche - position_rouge) / norm(position_blanche - position_rouge);
          r_ap = -r * n;
          r_bp = r * n;
          
          vr_avant = vitesse_avant_collision_blanche - vitesse_angulaire_avant_collision_rouge;
          J = calculerImpulsion(vr_avant, n, r_ap, r_bp);

          vitesse_apres_collision_blanche = vitesse_avant_collision_blanche + (J / m);
          vitesse_angulaire_apres_collision_blanche = vitesse_angulaire_avant_collision_blanche + inv(I) * cross(r_ap, J);
          
          vitesse_apres_collision_rouge = vitesse_avant_collision_rouge - (J / m);
          vitesse_angulaire_apres_collision_rouge = vitesse_angulaire_avant_collision_rouge - (inv(I) * cross(r_bp, J));
          
          q0_blanche = updateQ0ApresCollision(q0_blanche, vitesse_apres_collision_blanche, vitesse_angulaire_apres_collision_blanche);
          q0_rouge = updateQ0ApresCollision(q0_rouge, vitesse_apres_collision_rouge, vitesse_angulaire_apres_collision_rouge);

          if ((abs(norm(vitesse_apres_collision_blanche - cross(-vitesse_angulaire_apres_collision_blanche, r_pc)))) > 0.0005)
              % blanche_glisse = true;
          end
          if ((abs(norm(vitesse_apres_collision_rouge - cross(-vitesse_angulaire_apres_collision_rouge, r_pc)))) > 0.0005)
              % rouge_glisse = true;
          end
          
          col_rouge = 1;
        elseif (norm(position_blanche - position_jaune) <= r + r && ~block_blanche_jaune)
          blockers_blanc = resetBlockers();
          blockers_jaune = resetBlockers();
          block_blanche_jaune = true;
          n = (position_blanche - position_jaune) / norm(position_blanche - position_jaune);
          r_ap = -r * n;
          r_bp = r * n;
          
          vr_avant = vitesse_avant_collision_blanche - vitesse_angulaire_avant_collision_jaune;
          J = calculerImpulsion(vr_avant, n, r_ap, r_bp);

          vitesse_apres_collision_blanche = vitesse_avant_collision_blanche + (J / m);
          vitesse_angulaire_apres_collision_blanche = vitesse_angulaire_avant_collision_blanche + (inv(I) * cross(r_ap, J));
          
          vitesse_apres_collision_jaune = vitesse_avant_collision_jaune - (J / m);
          vitesse_angulaire_apres_collision_jaune = vitesse_angulaire_avant_collision_jaune - (inv(I) * cross(r_bp, J));
          
          q0_blanche = updateQ0ApresCollision(q0_blanche, vitesse_apres_collision_blanche, vitesse_angulaire_apres_collision_blanche);
          q0_jaune = updateQ0ApresCollision(q0_jaune, vitesse_apres_collision_jaune, vitesse_angulaire_apres_collision_jaune);
          
          if ((abs(norm(vitesse_apres_collision_blanche - cross(-vitesse_angulaire_apres_collision_blanche, r_pc)))) > 0.0005)
              % blanche_glisse = true;
          end
          if ((abs(norm(vitesse_apres_collision_jaune - cross(-vitesse_angulaire_apres_collision_jaune, r_pc)))) > 0.0005)
              % jaune_glisse = true;
          end

          col_jaune = 1;
        elseif (norm(position_rouge - position_jaune) <= r + r && ~block_jaune_rouge)
          blockers_rouge = resetBlockers();
          blockers_jaune = resetBlockers();
          block_jaune_rouge = true;
          
          n = (position_rouge - position_jaune) / norm(position_rouge - position_jaune);
          r_ap = -r * n;
          r_bp = r * n;
          
          vr_avant = vitesse_avant_collision_rouge - vitesse_angulaire_avant_collision_jaune;
          J = calculerImpulsion(vr_avant, n, r_ap, r_bp);

          vitesse_apres_collision_rouge = position_rouge + (J / m);
          vitesse_angulaire_apres_collision_rouge = vitesse_angulaire_avant_collision_rouge + (inv(I) * cross(r_ap, J));
          
          vitesse_apres_collision_jaune = vitesse_avant_collision_jaune - (J / m);
          vitesse_angulaire_apres_collision_jaune = vitesse_angulaire_avant_collision_jaune - (inv(I) * cross(r_bp, J));
          
          q0_rouge = updateQ0ApresCollision(q0_rouge, vitesse_apres_collision_rouge, vitesse_angulaire_apres_collision_rouge);
          q0_jaune = updateQ0ApresCollision(q0_jaune, vitesse_apres_collision_jaune, vitesse_angulaire_apres_collision_jaune);
          
          if ((abs(norm(vitesse_apres_collision_rouge - cross(-vitesse_angulaire_apres_collision_rouge, r_pc)))) > 0.0005)
              % rouge_glisse = true;
          end
          if ((abs(norm(vitesse_apres_collision_jaune - cross(-vitesse_angulaire_apres_collision_jaune, r_pc)))) > 0.0005)
              % jaune_glisse = true;
          end
          
        end
    end
    
    function [q0] = updateQ0ApresCollision(q0_intial, vitesse_apres_collision, vitesse_angulaire_apres_collision)
        q0 = q0_intial;
        q0(1) = vitesse_apres_collision(1);
        q0(2) = vitesse_apres_collision(2);
        q0(3) = vitesse_apres_collision(3);
        q0(7) = vitesse_angulaire_apres_collision(1);
        q0(8) = vitesse_angulaire_apres_collision(2);
        q0(9) = vitesse_angulaire_apres_collision(3);
    end
    
    function [blockers] = resetBlockers()
        % top, left, down, right
        blockers = [false, false, false, false];
        block_blanche_jaune = false;
        block_blanche_rouge = false;
        block_jaune_rouge = false;
    end
    % vitesse due au roulement: vr = -w x r_pc
    % tsi vitesse centre de masse > vitesse roulemen, objet glisse
    % vitesse glissement = vitesse - vr
    % Calcul methode de Runge-Kutta
    
    q0_blanche = [Vb0(1); Vb0(2); 0; xyb(1); xyb(2); r; 0; 0; 0];
    q0_jaune   = [0;    0;        0;      xyj(1); xyj(2); r; 0; 0; 0];
    q0_rouge   = [0;    0;        0;      xyr(1); xyr(2); r; 0; 0; 0];
    
    qs_blanche = q0_blanche;
    qs_jaune = q0_jaune;
    qs_rouge = q0_rouge;

    t0 = 0;
    
    % top, left, down, right
    blockers_blanc = resetBlockers();
    blockers_rouge = resetBlockers();
    blockers_jaune = resetBlockers();
    
    block_blanche_jaune = false;
    block_blanche_rouge = false;
    block_jaune_rouge = false;

    blanche_glisse = true;
    jaune_glisse = true;
    rouge_glisse = true;
    
    v_cm_blanche = [Vb0(1); Vb0(2); 0;];

    
    % Initialisation des variables de sorties
    tr = [t0];
    t  = [t0];
    xb = [q0_blanche(4)];
    yb = [q0_blanche(5)];
    xr = [q0_rouge(4)];
    yr = [q0_rouge(5)];
    xj = [q0_jaune(4)];
    yj = [q0_jaune(5)];

    while norm(v_cm_blanche) > 0.001
        
        q0_blanche = qs_blanche;
        q0_jaune   = qs_jaune;
        q0_rouge   = qs_rouge;
        
        position_blanche = [q0_blanche(4); q0_blanche(5); q0_blanche(6)];
        position_jaune   = [q0_jaune(4);   q0_jaune(5);   q0_jaune(6)];
        position_rouge   = [q0_rouge(4);   q0_rouge(5);   q0_rouge(6)];
        
        % check collisions
        % if collision set les nouvelles valeurs de q0(1); q0(2); q0(3); q0(5); q0(6); q0(7)
        % if collision reset glisse a true

        [q0_blanche, blockers_blanc, blanche_glisse] = checkCollisions(q0_blanche, blockers_blanc, blanche_glisse);
        [q0_rouge, blockers_rouge, rouge_glisse] = checkCollisions(q0_rouge, blockers_rouge, rouge_glisse);
        [q0_jaune, blockers_jaune, jaune_glisse] = checkCollisions(q0_jaune, blockers_jaune, jaune_glisse);

        [q0_blanche, q0_rouge, q0_jaune] = checkCollisionsBille(q0_blanche, q0_rouge, q0_jaune);
        
        v_cm_blanche         = [q0_blanche(1); q0_blanche(2); q0_blanche(3)];
        v_angulaire_blanche  = [q0_blanche(7); q0_blanche(8); q0_blanche(9)];
        
        v_cm_jaune           = [q0_jaune(1); q0_jaune(2); q0_jaune(3)];
        v_angulaire_jaune    = [q0_jaune(7); q0_jaune(8); q0_jaune(9)];
        
        v_cm_rouge           = [q0_rouge(1); q0_rouge(2); q0_rouge(3)];
        v_angulaire_rouge    = [q0_rouge(7); q0_rouge(8); q0_rouge(9)];
        
        v_roulement_blanche = cross(-v_angulaire_blanche, r_pc);
        v_roulement_jaune   = cross(-v_angulaire_jaune, r_pc);
        v_roulement_rouge 	= cross(-v_angulaire_rouge, r_pc);
        
        % determiner s'il ya glissement des boules
        if(abs(norm(v_cm_blanche - v_roulement_blanche)) <= 0.0005 && blanche_glisse)
            blanche_glisse = false;
            tr = [tr; t0];
        end
        
        if(abs(norm(v_cm_jaune) - norm(v_roulement_jaune)) <= 0.0005 && jaune_glisse)
            jaune_glisse = false;
        end
        
        if(abs(norm(v_cm_rouge) - norm(v_roulement_rouge)) <= 0.0005 && rouge_glisse)
            rouge_glisse = false;
        end
        
        qs_blanche = RK(q0_blanche, blanche_glisse, v_roulement_blanche, v_angulaire_blanche);
        qs_jaune = RK(q0_jaune, jaune_glisse, v_roulement_jaune, v_angulaire_jaune);
        qs_rouge = RK(q0_rouge, rouge_glisse, v_roulement_rouge, v_angulaire_rouge);
        
        t0 = t0 + dt;
        
        
        t = [t; t0];
        xb = [xb; q0_blanche(4)];
        yb = [yb; q0_blanche(5)];
        xr = [xr; q0_rouge(4)];
        yr = [yr; q0_rouge(5)];
        xj = [xj; q0_jaune(4)];
        yj = [yj; q0_jaune(5)];
    end
    
    temps_mis = t0;
    tr

    xb(end)
    yb(end)
    coll = col_rouge + col_jaune;
end