function RMatrixOZ = calc_rot_matrix_oz(angleRadians)
    
    cosAngleDimless = cos(angleRadians);
    sinAngleDimless = sin(angleRadians);
    
    RMatrixOZ = [cosAngleDimless,     -sinAngleDimless,   0.0;
                    sinAngleDimless,      cosAngleDimless,   0.0;
                         0.0,   0.0, 1.0];
end
