-- mesh sizes
num_cells = {25}

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

  -- 0 -- 1 -- 2
  coordinates = {
    {  0}, -- 0
    {0.5}, -- 1
    {  1}, -- 2
  },
      
  --------------------------------------
  -- the boundaries
  --------------------------------------
  boundaries = {
    {
      name = "left",
      connect = { {0} },
    },
    {
      name = "right",
      connect = { {2} },
      
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
      dimensions = {num_cells[1]},
      spacing = "linear",
      connect = {0, 1},
    },
    
    -- RIGHT BLOCK
    {
      dimensions = {num_cells[1]},
      spacing = "linear",
      --connect = {1, 2},
      connect = {2, 1},
    }
  } -- blocks
}

