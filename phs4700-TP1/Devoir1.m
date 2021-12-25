% Ismail Bakkouri - 1954157
% Mohammed Ariful Islam - 1950221
% John Maliha - 1984959
% Dawut Esse - 1956802

function [pcm MI aa] = 	Devoir1(pos, ar, va, Force)

matrice_rotation_x = [ 1, 0,       0;
                       0, cos(ar), -sin(ar);
                       0, sin(ar), cos(ar); ];  
                   
rayon = 1.5;
masse_vol = 4000;

% cylindre
longueur_cyl = 10;
x_cm_cyl = 0;
y_cm_cyl= 0;
z_cm_cyl= longueur_cyl / 2;
volume_cyl = longueur_cyl * pi * rayon^2;
masse_cyl = masse_vol * volume_cyl;

% cone
hauteur_cone = 3;
x_cm_cone = 0;
y_cm_cone = 0;
z_cm_cone = (hauteur_cone) / 4 + longueur_cyl;
volume_cone = (hauteur_cone * pi * rayon^2) / 3;
masse_cone = masse_vol * volume_cone;

% aileron
epaisseur_aileron = 0.15;
hauteur_aileron = 1.5;
largeur_aileron = 2.5;
masse_aileron = 2000;

x_cm_aileron1 = -largeur_aileron / 2 - rayon;
x_cm_aileron2 = 0;
x_cm_aileron3 = largeur_aileron / 2 + rayon;
x_cm_aileron4 = 0;

y_cm_aileron1 = 0;
y_cm_aileron2 = largeur_aileron / 2 + rayon;
y_cm_aileron3 = 0;
y_cm_aileron4 = -largeur_aileron / 2 - rayon;

z_cm_aileron = hauteur_aileron / 2;

masse_total_systeme = masse_cyl + masse_cone + masse_aileron * 4;

% composantes du centre de masse
x_cm = (x_cm_cyl * masse_cyl + x_cm_cone * masse_cone + x_cm_aileron1 * masse_aileron + x_cm_aileron2 * masse_aileron + x_cm_aileron3 * masse_aileron + x_cm_aileron4 * masse_aileron) / masse_total_systeme;
y_cm = (y_cm_cyl * masse_cyl + y_cm_cone * masse_cone + y_cm_aileron1 * masse_aileron + y_cm_aileron2 * masse_aileron + y_cm_aileron3 * masse_aileron + y_cm_aileron4 * masse_aileron) / masse_total_systeme;
z_cm = (z_cm_cyl * masse_cyl + z_cm_cone * masse_cone + (z_cm_aileron * masse_aileron * 4)) / masse_total_systeme;


% PCM
pcm_avant_rotation_translation = transpose([x_cm, y_cm, z_cm]);
pcm = (matrice_rotation_x * pcm_avant_rotation_translation) + pos;


% Moment d'inertie cylindre par rapport au centre de masse du cylindre
matrice_inertie_cylindre = [masse_cyl / 4 * rayon^2 + (masse_cyl / 12) * longueur_cyl^2, 0, 0;
         0, masse_cyl / 4 * rayon^2 + masse_cyl / 12 * longueur_cyl^2, 0;
         0, 0, masse_cyl / 2 * rayon^2 ];
     
% Moment d'inertie cone par rapport au centre de masse du cone
matrice_inertie_cone = [ masse_cone*(12 * rayon^2 + 3 * hauteur_cone^2)/80, 0, 0;
         0, masse_cone * (12 * rayon^2 + 3 * hauteur_cone^2) / 80, 0;
         0, 0, masse_cone * 3 * rayon^2 / 10 ];
     
% Moment d'inertie ailerons par rapport au cm des ailerons
    function [Matrice_inertie_finale] = MI_aileron(a, b ,c, m, cm_aileron)
        Matrice_inertie = m*[(b^2+c^2)/12, 0, 0;
                          0,(a^2+c^2)/12, 0;
                          0, 0,(a^2+b^2)/12;];
        Matrice_inertie_finale = translater_MI(m, cm_aileron, Matrice_inertie);
                            
    end

    function [Matrice_inertie_translatee] = translater_MI(m, cm_obj, MI_obj)
        d = pcm_avant_rotation_translation - cm_obj;
        matrice_T = [ d(2)^2 + d(3)^2, -1 * d(1) * d(2), -1 * d(1) * d(3);
                      -1 * d(2) * d(1), d(1)^2 + d(3)^2, -1 * d(2) * d(3);
                      -1 * d(3) * d(1), -1 * d(3) * d(2), d(1)^2 + d(2)^2; ];
        
        Matrice_inertie_translatee = MI_obj + (m * matrice_T);
    end

MI_aileron_1_translatee = MI_aileron(0.15, 2.5, 1.5, masse_aileron,[x_cm_aileron1; y_cm_aileron1; z_cm_aileron]);
MI_aileron_2_translatee = MI_aileron(2.5, 0.15, 1.5, masse_aileron,[x_cm_aileron2; y_cm_aileron2; z_cm_aileron]);
MI_aileron_3_translatee = MI_aileron(0.15, 2.5, 1.5, masse_aileron,[x_cm_aileron3; y_cm_aileron3; z_cm_aileron]);
MI_aileron_4_translatee = MI_aileron(2.5, 0.15, 1.5, masse_aileron,[x_cm_aileron4; y_cm_aileron4; z_cm_aileron]);

MI_cylindre_translatee = translater_MI(masse_cyl,[x_cm_cyl; y_cm_cyl; z_cm_cyl], matrice_inertie_cylindre);
MI_cone_translatee = translater_MI(masse_cone,[x_cm_cone; y_cm_cone; z_cm_cone], matrice_inertie_cone);

MI_translatee = MI_aileron_1_translatee + MI_aileron_2_translatee + MI_aileron_3_translatee + MI_aileron_4_translatee + MI_cylindre_translatee + MI_cone_translatee;

MI = matrice_rotation_x * MI_translatee * transpose(matrice_rotation_x);

% Moment de force
tetha = Force(2);
phi = Force(3);
force_decomposee = Force(1) * [cos(phi) * sin(tetha); sin(phi) * sin(tetha); cos(tetha)];
Torque = cross((pos - pcm), force_decomposee);

% Moment cinetique
L =  MI * va;

% Acceleration angulaire

aa = MI\(Torque + cross(L, va));

end