function [xi yi zi couleur] = Devoir4(Robs, xyp, he)
    Rp = 0.005;       % rayon paille
    Hp = 0.25;        % hauteur paille
    Rv = 0.04;        % rayon verre
    Hv = 0.2;         % hauteur verre
    n_eau = 1.333;    % indice de refraction eau
    n_air = 1;        % indice de refraction verre

    angleCritique = abs(asin(n_eau/n_air));

    anglePolaireMax = 359;
    angleAzimutalMax = 359;
    r_paille = [xyp(1); xyp(2); Hp/2];
    hauteur_section = Hp / 5;
    paille = Paille(r_paille, hauteur_section, Rp);
    eau = Cylindre([0; 0; he / 2], he, Rv);
    
    couvercle = Disque([0; 0; he], [0; 0; 1], Rv);
    base = Disque([0; 0; 0], [0;0;-1], Rv);
    
    couvercle_paille = Disque([xyp(1); xyp(2); Hp], [0; 0; 1], Rp);
    base_paille = Disque([xyp(1); xyp(2); 0], [0; 0; -1], Rp);
  
    function rayonResultant = calculerNouveauRayon(rayon, normale, point_intersection, n1, n2) 
      j = cross(rayon.direction, normale)/norm(cross(rayon.direction, normale));
      k = cross(normale, j)/norm(cross(normale, j));
      si = dot(rayon.direction, k);
      angleIncident = asin(si);

      if((angleIncident > angleCritique || angleIncident < -angleCritique))
        % Reflexion
        direction = normale * sqrt(1 - (si * si)) + (k * si);
      else
        % Refraction
        interieur = ~interieur;
        st = (n1/n2) * si;
        direction = -normale * sqrt(1 - (st * st)) + k * st;
      end
      
      direction = direction/norm(direction);
      rayonResultant = Rayon(point_intersection, direction);
    end

    function premiereCollision = trouverPremiereIntersection(distances)
        index = -1;
        dist = Inf;
        for i = 1:length(distances)
            if(dist >= distances(i) && distances(i) > 0.0001)
                index = i;
                dist = distances(i);
            end
        end
        premiereCollision = index;
    end

  xi = [];
  yi = [];
  zi = [];
  couleur = [];
  
  for anglePolaire = 0:0.3:anglePolaireMax
    for angleAzimutal = 0:0.3:angleAzimutalMax
        
      direction_initiale = [sin(anglePolaire) * cos(angleAzimutal); sin(anglePolaire) * sin(angleAzimutal); cos(anglePolaire)];
      direction_initiale = direction_initiale / norm(direction_initiale);
      rayon_initial = Rayon(Robs, direction_initiale);
      rayon_courant = rayon_initial;
      nb_bounces = 0;
      distance_totale = 0;
      interieur = false;
        
      while(nb_bounces < 10)
          distances = [];
          intersection = [];

          [intersect, point_intersection, distance] = eau.intersection(rayon_courant);
          distances = [distances; distance];
          intersection = [intersection; point_intersection];
         
          [intersect, point_intersection, distance] = couvercle.intersection(rayon_courant);
          distances = [distances; distance];
          intersection = [intersection; point_intersection];
      
          [intersect, point_intersection, distance] = base.intersection(rayon_courant);
          distances = [distances; distance];
          intersection = [intersection; point_intersection];   
          
          [intersect, point_intersection, distance] = paille.section_rouge.intersection(rayon_courant);
          distances = [distances; distance];
          intersection = [intersection; point_intersection];
 
          [intersect, point_intersection, distance] = paille.section_orange.intersection(rayon_courant);
          distances = [distances; distance];
          intersection = [intersection; point_intersection];
 
          [intersect, point_intersection, distance] = paille.section_magenta.intersection(rayon_courant);
          distances = [distances; distance];
          intersection = [intersection; point_intersection];
 
          [intersect, point_intersection, distance] = paille.section_vert.intersection(rayon_courant);
          distances = [distances; distance];
          intersection = [intersection; point_intersection];
 
          [intersect, point_intersection, distance] = paille.section_bleu.intersection(rayon_courant);
          distances = [distances; distance];
          intersection = [intersection; point_intersection];
          
          [intersect, point_intersection, distance] = couvercle_paille.intersection(rayon_courant);
          distances = [distances; distance];
          intersection = [intersection; point_intersection];
          
          [intersect, point_intersection, distance] = base_paille.intersection(rayon_courant);
          distances = [distances; distance];
          intersection = [intersection; point_intersection];

          % 1  = cylindre eau
          % 2  = couvercle eau
          % 3  = base eayu
          % 4  = section rouge
          % 5  = section orange
          % 6  = section magenta
          % 7  = section vert
          % 8  = section bleu
          % 9  = couvercle paille
          % 10 = base paille
          
          indexDistance = trouverPremiereIntersection(distances);
          if(indexDistance < 0)
              % aucune intersection, continuer l'iteration
            break;
          end
          
          point = [intersection(indexDistance * 3 - 2); intersection(indexDistance * 3 - 1); intersection(indexDistance * 3)];

          if(indexDistance < 4)
            distance_totale = distance_totale + distances(indexDistance);
            
            if(indexDistance == 1) % cote du cylindre
                normale = eau.getNormale(point, interieur);
            elseif(indexDistance == 2) % couvercle du cylindre
                normale = couvercle.getNormale(interieur);
            elseif(indexDistance == 3) % base du cylindre
                normale = base.getNormale(interieur);
            end
              
            if(interieur)
                 nouveau_rayon = calculerNouveauRayon(rayon_courant, normale, point, n_eau, n_air);
                 rayon_courant = nouveau_rayon;
            else
                 nouveau_rayon = calculerNouveauRayon(rayon_courant, normale, point, n_air, n_eau);
                 rayon_courant = nouveau_rayon;
            end
          
            nb_bounces = nb_bounces + 1;
          else
              % touche la paille
              distance_totale = distance_totale + distances(indexDistance);
              point_final = rayon_initial.trouverPoint(distance_totale);

              xi = [xi;  point_final(1)];
              yi = [yi;  point_final(2)];
              zi = [zi;  point_final(3)];
              if(indexDistance == 4 || indexDistance == 10) % rouge & base
                  couleur = [couleur; 1];
              elseif(indexDistance == 5) % orange
                  couleur = [couleur; 2];
              elseif(indexDistance == 6) % magenta
                  couleur = [couleur; 3];
              elseif(indexDistance == 7) % vert
                  couleur = [couleur; 4];
              elseif(indexDistance == 8 || indexDistance == 9) % bleu & couvercle
                  couleur = [couleur; 5];
              end
            break;
          end
      end
    end
    anglePolaire
  end
end