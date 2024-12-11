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
  type = "auto",
  
  -- the boundaries
  boundaries = {
    {
      name = "left",
      id = 0,
    },
    {
      name = "right",
      id = 1,
      
      where = function(x)
        return true
      end
    },
  }
}

