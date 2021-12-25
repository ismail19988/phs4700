classdef Disque
  properties
    normale = [];
    centre = [];
    rayon = 0;
  end

  methods
    function obj = Disque(centre, normale, rayon)
       obj.normale = normale;
       obj.centre = centre;
       obj.rayon = rayon;
    end

    function normale =  getNormale(obj, interieur)
      if (interieur)
        normale = -obj.normale;
      else
        normale = obj.normale;
      end
    end

function [intersect, point_intersection, distance] = intersection(obj, rayon)
      intersect = false;
      point_intersection = [0; 0; 0];
      distance = 0;
      
      if(rayon.point_depart(3) > obj.centre(3))
        normaleCalcul = [0; 0; -1];
      else
        normaleCalcul = [0; 0; 1];
      end
      
      A = dot(normaleCalcul, rayon.direction);
      B = dot(normaleCalcul, rayon.point_depart - obj.centre);

      if (A > 0.0001)
          t = -B/A;
          if (t > 0)
            point = rayon.trouverPoint(t) - obj.centre;
            if (norm(point) <= obj.rayon)
              intersect = true;
              point_intersection = rayon.trouverPoint(t);
              distance = norm(point_intersection - rayon.point_depart);
            end
          end
      end
    end
  end

end