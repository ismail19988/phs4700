classdef Cylindre
  properties
    r_cylindre = [];
    r_a1 = [];
    r_a2 = [];
    s = [];
    hauteur = 0;
    rayon = 0;
    cm = []
  end
  
  methods
    function obj = Cylindre(r_cylindre, hauteur, rayon)
      obj.r_cylindre = r_cylindre;
      obj.hauteur = hauteur;
      obj.rayon = rayon;
      obj.r_a1 = r_cylindre;
      obj.r_a1(3) = obj.r_a1(3) - hauteur / 2;
      
      obj.r_a2 = r_cylindre;
      obj.r_a2(3) = obj.r_a2(3) + hauteur / 2;
      
      obj.s = (obj.r_a2 - obj.r_a1) / norm(obj.r_a2 - obj.r_a1);
      obj.cm = [0; 0; hauteur/2];
    end
    
    function [intersect, point_intersection, distance] = intersection(obj, rayon)
      point_intersection = [0; 0; 0];
      intersect = false;
      distance = 0;
      r_a0 = cross(cross(obj.s, rayon.point_depart - obj.r_a1), obj.s);
      v_a = cross(cross(obj.s, rayon.direction), obj.s);
      
      A = dot(v_a, v_a);
      B = 2 * dot(r_a0, v_a);
      C = dot(r_a0, r_a0) - obj.rayon * obj.rayon;
      
      [valide, t] = obj.evaluerFormuleQuadratic(A, B, C);
      
      if (valide)
        point_intersection = rayon.trouverPoint(t);
        EntreCylindre = obj.estEntreIntersection(point_intersection);
        if (EntreCylindre && t > 0)
          intersect = true;
          distance = norm(point_intersection - rayon.point_depart);
        end
      end
    end
    
    function entreIntersection = estEntreIntersection(obj, r_intersection)
      interieur_r_a1 = dot(r_intersection - obj.r_a1, obj.s);
      interieur_r_a2 = dot(r_intersection - obj.r_a2, obj.s);
      entreIntersection = interieur_r_a1 > 0 && interieur_r_a2 < 0;
    end
    
    function normale = getNormale(obj, point, interieur)
        normale = point - obj.cm;
        normale(3) = 0;
        normale = normale / norm(normale);
        if(interieur)
            % si nous somme dans le cylindre, la normal doit etre dans le
            % bon sens (pointer vers l'intérieur)
            normale = -1 * normale;
        end
    end
    
    function [valide, t] = evaluerFormuleQuadratic(obj, A, B, C)
      if (A < 0.0001)
        if (B < 0.0001)
          valide = false;
          t = 0;
        else
          valide = true;
          t = -B/A;
          if (t < 0)
            valide = false;
          end
        end
      else 
        discriminant = B*B - 4*A*C;
        if (discriminant < 0)
          valide = false;
          t = 0;
        else
          valide = true;
          t = 0;
          t_1 = (-B - sqrt(discriminant)) / (2*A);
          t_2 = (-B + sqrt(discriminant)) / (2*A);
          if (t_1 < 0 && t_2 < 0)
            valide = false;
          elseif (t_1 > 0 && t_2 < 0)
            t = t_1;
          elseif (t_1 < 0 && t_2 > 0)
            t = t_2;
          elseif (t_1 > 0 && t_2 > 0)
            if (t_1 < t_2)
              t = t_1;
            else 
              t = t_2;
            end
          end
        end
      end
    end
    
  end
end