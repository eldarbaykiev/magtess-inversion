function RMatrixOY = calc_rot_matrix_oy(angleRadians)
    
    cosAngleDimless = cos(angleRadians);
    sinAngleDimless = sin(angleRadians);
    
    RMatrixOY = [   cosAngleDimless      0.0   sinAngleDimless;
                    0.0                  1.0   0.0;
                    -sinAngleDimless     0.0   cosAngleDimless];
end

