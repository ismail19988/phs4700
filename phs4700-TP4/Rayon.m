classdef Rayon
   properties
     point_depart = [];
     direction = [];
   end
   
   methods
     function obj = Rayon(point_depart, direction)
       obj.point_depart = point_depart;
       obj.direction = direction;
     end
     
     function point = trouverPoint(obj, t)
         point =  obj.point_depart + (t * obj.direction);
     end
      
     function distance = trouverDistance(obj, point)
         distance = norm(point - obj.point_depart);
     end
   end
end