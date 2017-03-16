function vect2 = aux_vect_rot(vect1, lon1, lat1, lon2, lat2)
    b1 = conv_lon2sph(lon1);
    a1 = conv_lat2sph(lat1);
    b2 = conv_lon2sph(lon2);
    a2 = conv_lat2sph(lat2);
    
    E = [-1 0 0;
          0 1 0;
          0 0 1]; %NORTH EAST UP
    
    vect2 = (calc_rot_matrix_oz(b2) * calc_rot_matrix_oy(pi/2.0-a2) * E) \ calc_rot_matrix_oz(b1) * calc_rot_matrix_oy(pi/2.0-a1) * E * vect1;
    return;
end