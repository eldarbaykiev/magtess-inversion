function RMatrixOX = calc_rot_matrix_ox(angleRadians)
    
    cosAngleDimless = cos(angleRadians);
    sinAngleDimless = sin(angleRadians);
    
    RMatrixOX = [   1.0,       0.0,                   0.0;
                    0.0,       cosAngleDimless,     -sinAngleDimless;
                    0.0,       sinAngleDimless,      cosAngleDimless];
end