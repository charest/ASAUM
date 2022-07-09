-- mesh sizes
num_cells = {25, 10, 10}

--------------------------------------------------------------------------------
-- Begin Main Input
--------------------------------------------------------------------------------

driver = {
  -- The case prefix and postfixes
  output_prefix = "sod",
  output_formats = "csv,vtk",
}

vtk = {
  encoding = "ascii"
}

--------------------------------------------------------------------------------
-- Mesh Input
--------------------------------------------------------------------------------

mesh = {
  type = "structured",

  -- 5 -- 4 -- 3
  -- |    |    |
  -- 0 -- 1 -- 2
  coordinates = {
    {  0, -0.05, -0.05}, -- 0
    {0.5, -0.05, -0.05}, -- 1
    {  1, -0.05, -0.05}, -- 2
    {  1,  0.05, -0.05}, -- 3
    {0.5,  0.05, -0.05}, -- 4
    {  0,  0.05, -0.05}, -- 5
    {  0, -0.05,  0.05}, -- 6
    {0.5, -0.05,  0.05}, -- 7
    {  1, -0.05,  0.05}, -- 8
    {  1,  0.05,  0.05}, -- 9
    {0.5,  0.05,  0.05}, -- 10
    {  0,  0.05,  0.05}, -- 11
  },

  --rotate = {-1, 3, 2},

  -- the boundaries
  boundaries = {
    {
      name = "left",
      connect = {{0, 5, 11, 6}},
      type = "extrapolate",
    },
    {
      name = "right",
      connect = {{2, 8, 9, 3}},
      type = "extrapolate",
      
      where = function(x)
        return true
      end
    },
    {
      name = "bottom",
      connect = {{0, 6, 7, 1}, {1, 7, 8, 2}},
      type = "symmetry",
    },
    {
      name = "top",
      connect = {{5, 4, 10, 11}, {4, 3, 9, 10}},
      type = "symmetry",
      
      where = function(x)
        return true
      end
    },
    {
      name = "south",
      connect = {{0, 1, 4, 5}, {1, 2, 3, 4}},
      type = "symmetry",
    },
    {
      name = "north",
      connect = {{11, 10, 7, 6}, {10, 9, 8, 7}},
      type = "symmetry",
      
      where = function(x)
        return true
      end
    }
  }, -- boundaries


  --------------------------------------
  -- The Blocks
  --------------------------------------
  blocks = {
    -- LEFT BLOCK
    {
      dimensions = {num_cells[1], num_cells[2], num_cells[3]},
      spacing = "linear",
      connect = {0, 1, 4, 5, 6, 7, 10, 11},
    },
    
    -- RIGHT BLOCK
    {
      dimensions = {num_cells[1], num_cells[2], num_cells[3]},
      spacing = "linear",
      -- aligned - x
      --connect = {1, 2, 3, 4, 7, 8, 9, 10},
      -- rotate once - x
      --connect = {4, 3, 9, 10, 1, 2, 8, 7},
      -- rotate twice - x
      --connect = {10, 9, 8, 7, 4, 3, 2, 1},
      -- rotate thrice - x
      connect = {7, 8, 2, 1, 10, 9, 3, 4},
    }
  } -- blocks
}

