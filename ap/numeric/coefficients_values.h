/* GSL-style data file */
#ifdef USE_ROW_MAJOR

        double coefficients[21*21] = {
                1, 0, 0, 0, 0, 0, -10, 0, 0, -10, 15, 0, -30, 0, 15, -6, 0, 30,
                30, 0, -6,
                0, 0, 0, 0, 0, 0, 10, 0, 0, 0, -15, 0, 15, 0, 0, 6, 0, -15,
                -15, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 10, 0, 0, 15, 0, -15, 0, 0, -15,
                -15, 0, 6,
                0, 1, 0, 0, 0, 0, -6, 0, -11, 0, 8, 0, 10, 18, 0, -3, 0, 1,
                -10, -8, 0,
                0, 0, 1, 0, 0, 0, 0, -11, 0, -6, 0, 18, 10, 0, 8, 0, -8, -10,
                1, 0, -3,
                0, 0, 0, 0, 0, 0, -4, 0, 0, 0, 7, 0, -3.5, 0, 0, -3, 0, 3.5,
                3.5, 0, 0,
                0, 0, 0, 0, 0, 0, 0, -5, 0, 0, 0, 14, 18.5, 0, 0, 0, -8, -18.5,
                -13.5, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, -5, 0, 0, 0, 18.5, 14, 0, 0, 0, -13.5,
                -18.5, -8, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 0, 0, -3.5, 0, 7, 0, 0, 3.5,
                3.5, 0, -3,
                0, 0, 0, 0.5, 0, 0, -1.5, 0, 0, 0, 1.5, 0, -1.5, 0, 0, -0.5, 0,
                1.5, 1, 0, 0,
                0, 0, 0, 0, 1, 0, 0, -4, -4, 0, 0, 5, 10, 5, 0, 0, -2, -6, -6,
                -2, 0,
                0, 0, 0, 0, 0, 0.5, 0, 0, 0, -1.5, 0, 0, -1.5, 0, 1.5, 0, 0, 1,
                1.5, 0, -0.5,
                0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, -1, 0, 0.25, 0, 0, 0.5, 0, -0.25,
                -0.25, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -3, -3.5, 0, 0, 0, 2, 3.5,
                2.5, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.25, 0, 0, 0, 0, -0.75,
                -1.25, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1.25, 0, 0, 0, 0, -1.25,
                -0.75, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -3.5, -3, 0, 0, 0, 2.5,
                3.5, 2, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0.25, 0, -1, 0, 0, -0.25,
                -0.25, 0, 0.5,
                0, 0, 0, 0, 0, 0, 0, 16, 0, 0, 0, -32, -32, 0, 0, 0, 16, 32,
                16, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, -16, 0, 0, 0, 32, 32, 0, 0, 0, -16,
                -32, -16, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8*SQRT2, 0, 0, 0, 0,
                -8*SQRT2, -8*SQRT2, 0, 0};
#else

        double coefficients[21*21] = {
                1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                -10, 10, 0, -6, 0, -4, 0, 0, 0, -1.5, 0, 0, 0.5, 0, 0, 0, 0, 0,
                0, 0, 0,
                0, 0, 0, 0, -11, 0, -5, 0, 0, 0, -4, 0, 0, 1, 0, 0, 0, 0, 16, 0,
                0,
                0, 0, 0, -11, 0, 0, 0, -5, 0, 0, -4, 0, 0, 0, 0, 0, 1, 0, 0, -16,
                0,
                -10, 0, 10, 0, -6, 0, 0, 0, -4, 0, 0, -1.5, 0, 0, 0, 0, 0, 0.5,
                0, 0, 0,
                15, -15, 0, 8, 0, 7, 0, 0, 0, 1.5, 0, 0, -1, 0, 0, 0, 0, 0, 0,
                0, 0,
                0, 0, 0, 0, 18, 0, 14, 0, 0, 0, 5, 0, 0, -3, 0, 0, 0, 0, -32, 0,
                0,
                -30, 15, 15, 10, 10, -3.5, 18.5, 18.5, -3.5, -1.5, 10, -1.5,
                0.25, -3.5, 1.25, 1.25, -3.5, 0.25, -32, 32, 8*SQRT2,
                0, 0, 0, 18, 0, 0, 0, 14, 0, 0, 5, 0, 0, 0, 0, 0, -3, 0, 0, 32,
                0,
                15, 0, -15, 0, 8, 0, 0, 0, 7, 0, 0, 1.5, 0, 0, 0, 0, 0, -1, 0, 0,
                0,
                -6, 6, 0, -3, 0, -3, 0, 0, 0, -0.5, 0, 0, 0.5, 0, 0, 0, 0, 0, 0,
                0, 0,
                0, 0, 0, 0, -8, 0, -8, 0, 0, 0, -2, 0, 0, 2, 0, 0, 0, 0, 16, 0,
                0,
                30, -15, -15, 1, -10, 3.5, -18.5, -13.5, 3.5, 1.5, -6, 1, -0.25,
                3.5, -0.75, -1.25, 2.5, -0.25, 32, -16, -8*SQRT2,
                30, -15, -15, -10, 1, 3.5, -13.5, -18.5, 3.5, 1, -6, 1.5, -0.25,
                2.5, -1.25, -0.75, 3.5, -0.25, 16, -32, -8*SQRT2,
                0, 0, 0, -8, 0, 0, 0, -8, 0, 0, -2, 0, 0, 0, 0, 0, 2, 0, 0, -16,
                0,
                -6, 0, 6, 0, -3, 0, 0, 0, -3, 0, 0, -0.5, 0, 0, 0, 0, 0, 0.5, 0,
                0, 0};
#endif
