# expand_matrix_list

    Code
      expand_matrix_list(list(A, B))
    Output
      15 x 15 sparse Matrix of class "dgCMatrix"
                                                      
       [1,] 1 5  9 13 17  .  .  .  .  .  .  .  .  .  .
       [2,] 2 6 10 14 18  .  .  .  .  .  .  .  .  .  .
       [3,] 3 7 11 15 19  .  .  .  .  .  .  .  .  .  .
       [4,] 4 8 12 16 20  .  .  .  .  .  .  .  .  .  .
       [5,] . .  .  .  . 21 25 29 33 37 41 45 49 53 57
       [6,] . .  .  .  . 22 26 30 34 38 42 46 50 54 58
       [7,] . .  .  .  . 23 27 31 35 39 43 47 51 55 59
       [8,] . .  .  .  . 24 28 32 36 40 44 48 52 56 60
       [9,] . .  .  .  .  .  .  .  .  .  .  .  .  .  .
      [10,] . .  .  .  .  .  .  .  .  .  .  .  .  .  .
      [11,] . .  .  .  .  .  .  .  .  .  .  .  .  .  .
      [12,] . .  .  .  .  .  .  .  .  .  .  .  .  .  .
      [13,] . .  .  .  .  .  .  .  .  .  .  .  .  .  .
      [14,] . .  .  .  .  .  .  .  .  .  .  .  .  .  .
      [15,] . .  .  .  .  .  .  .  .  .  .  .  .  .  .

# expand_matrix

    Code
      expand_matrix(A)
    Output
      5 x 5 sparse Matrix of class "dgCMatrix"
                       
      [1,] 1 5  9 13 17
      [2,] 2 6 10 14 18
      [3,] 3 7 11 15 19
      [4,] 4 8 12 16 20
      [5,] . .  .  .  .

# expand_matrix_symmetric

    Code
      expand_matrix_symmetric(A)
    Output
      9 x 9 sparse Matrix of class "dgCMatrix"
                                    
       [1,]  .  .  .  . 1 5  9 13 17
       [2,]  .  .  .  . 2 6 10 14 18
       [3,]  .  .  .  . 3 7 11 15 19
       [4,]  .  .  .  . 4 8 12 16 20
       [5,]  1  2  3  4 . .  .  .  .
       [6,]  5  6  7  8 . .  .  .  .
       [7,]  9 10 11 12 . .  .  .  .
       [8,] 13 14 15 16 . .  .  .  .
       [9,] 17 18 19 20 . .  .  .  .

# trace_neighbors

    Code
      trace_neighbors(matrix(locs[10, ], 1, 2), grid, dat_neighbors)
    Output
      [1] 946

# get_neighbors

    Code
      dat_neighbors
    Output
      [[1]]
      # A tibble: 12,344 x 2
         grid_cell grid_cell_neighbors
             <int>               <int>
       1         1                   1
       2         1                   2
       3         1                   3
       4         1                  23
       5         1                  24
       6         1                  25
       7         1                  45
       8         1                  46
       9         1                  47
      10         2                   1
      # ... with 12,334 more rows
      
      [[2]]
      # A tibble: 11,998 x 2
         grid_cell grid_cell_neighbors
             <int>               <int>
       1         1                   1
       2         1                   2
       3         1                   3
       4         1                  43
       5         1                  44
       6         1                  45
       7         1                  85
       8         1                  86
       9         1                  87
      10         2                   1
      # ... with 11,988 more rows
      
      [[3]]
      # A tibble: 43,866 x 2
         grid_cell grid_cell_neighbors
             <int>               <int>
       1         1                   1
       2         1                   2
       3         1                   3
       4         1                  83
       5         1                  84
       6         1                  85
       7         1                 165
       8         1                 166
       9         1                 167
      10         2                   1
      # ... with 43,856 more rows
      

# get_neighbors_layers

    Code
      dat_neighbors
    Output
      [[1]]
      # A tibble: 3,036 x 3
         grid_cell_neighbors grid_cell layer
                       <int>     <int> <chr>
       1                   1         1 1~2  
       2                   2         1 1~2  
       3                  23         1 1~2  
       4                  24         1 1~2  
       5                   1         2 1~2  
       6                   2         2 1~2  
       7                   3         2 1~2  
       8                   4         2 1~2  
       9                  23         2 1~2  
      10                  24         2 1~2  
      # ... with 3,026 more rows
      
      [[2]]
      # A tibble: 11,156 x 3
         grid_cell_neighbors grid_cell layer
                       <int>     <int> <chr>
       1                   1         1 2~3  
       2                   2         1 2~3  
       3                  43         1 2~3  
       4                  44         1 2~3  
       5                   1         2 2~3  
       6                   2         2 2~3  
       7                   3         2 2~3  
       8                   4         2 2~3  
       9                  43         2 2~3  
      10                  44         2 2~3  
      # ... with 11,146 more rows
      

# calc_r_grid

    Code
      calc_r_grid(grid)
    Output
      [1] 1.8856181 0.8931875

