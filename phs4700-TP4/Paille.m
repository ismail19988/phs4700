classdef Paille
    properties
      section_rouge;
      section_orange;
      section_magenta;
      section_vert;
      section_bleu;
      sections = [];
    end
    
    methods
      function obj = Paille(r_paille, hauteur_section_paille, rayon_paille)
        obj.section_rouge   = Cylindre([r_paille(1); r_paille(2); hauteur_section_paille / 2],      hauteur_section_paille, rayon_paille);
        obj.section_orange  = Cylindre([r_paille(1); r_paille(2); 3 * hauteur_section_paille / 2],  hauteur_section_paille, rayon_paille);
        obj.section_magenta = Cylindre([r_paille(1); r_paille(2); 5 * hauteur_section_paille / 2],  hauteur_section_paille, rayon_paille);
        obj.section_vert    = Cylindre([r_paille(1); r_paille(2); 7 *hauteur_section_paille / 2],   hauteur_section_paille, rayon_paille);
        obj.section_bleu    = Cylindre([r_paille(1); r_paille(2); 9 * hauteur_section_paille / 2],  hauteur_section_paille, rayon_paille);
      end
    end
end