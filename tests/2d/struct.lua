-- mesh sizes
num_cells = {25, 10}

--------------------------------------------------------------------------------
-- Begin Main Input
--------------------------------------------------------------------------------

driver = {
  -- The case prefix and postfixes
  output_prefix = "out",
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
    {  0, -0.05}, -- 0
    {0.5, -0.05}, -- 1
    {  1, -0.05}, -- 2
    {  1,  0.05}, -- 3
    {0.5,  0.05}, -- 4
    {  0,  0.05}, -- 5
  },
      
  --------------------------------------
  -- the boundaries
  --------------------------------------
  boundaries = {
    {
      name = "left",
      connect = { {0,5} },
    },
    {
      name = "right",
      connect = { {2,3} },
      
      where = function(x)
        return true
      end
    },
    {
      name = "bottom",
      connect = { {0,1}, {1,2} },
    },
    {
      name = "top",
      connect = { {3,4}, {4,5} },
      
      where = function(x)
        return true
      end
    },
  }, -- boundaries


  --------------------------------------
  -- The Blocks
  --------------------------------------
  blocks = {
    -- LEFT BLOCK
    {
      dimensions = {num_cells[1], num_cells[2]},
      spacing = "linear",
      connect = {0, 1, 4, 5},
    },
    
    -- RIGHT BLOCK
    {
      dimensions = {num_cells[1], num_cells[2]},
      spacing = "linear",
      --connect = {1, 2, 3, 4},
      connect = {3, 4, 1, 2},
    }
  } -- blocks

}

