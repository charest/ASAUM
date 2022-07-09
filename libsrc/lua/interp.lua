--------------------------------------------------------------------------------
-- very similar to matlabs interp1 function
--------------------------------------------------------------------------------
function interp1(x, y, q, method)

  method = method or "linear"
 
  assert(#x == #y, "length of x and y must match")

  first = 1
  last = #y
  middle = math.floor( (first+last)/2 )

  --
  -- Binary search
  --
  while (first <= last) do
    if x[middle] < q then
      first = middle + 1
    elseif x[middle] == q then
      break
    else
      last = middle - 1
    end
    middle = math.floor( (first + last) / 2 )
  end

  --
  -- no exact match
  --
  if first > last then
  
    -- linear interpolation
    if method == "linear" then
      local res = {}
      local fact = (q - x[last]) / (x[first] - x[last])
      local mincols = math.min(#y[first], #y[last])
      for i = 1, mincols do 
        local val = y[last][i] + fact*(y[first][i]-y[last][i])
        table.insert(res, val)
      end 
      return res
    -- nearest neighbor solution
    elseif method == "nearest" then
      local left = math.abs(q - x[last])
      local right = math.abs(x[first] - q)
      if left < right then
        return y[last]
      else
        return y[first]
      end
    -- return the next value within the bracket
    elseif method == "next" then
      return y[first]
    -- return the previous value within the bracket
    elseif method == "previous" then
      return y[last]
    -- error
    else
      error("Uknown interp method "..method)
    end

  --
  -- an exact match
  --
  else

    return y[middle]

  end

end

--------------------------------------------------------------------------------
-- very similar to matlabs interp1 function
--------------------------------------------------------------------------------
function table_interp1(tbl, q, col, method)

  method = method or "linear"

  col = col or 1
 
  if col and col < 1 then
    col = nil
    error("Column must be greater than 0")
  end

  first = 1
  last = #tbl
  middle = math.floor( (first+last)/2 )

  local function exclude_col(t, c)
    -- result will be a table
    if #t > 2 then
      local res = {}
      for i,v in ipairs(t) do
        if i ~= c then
          table.insert(res, v)
        end
      end
      return res
    -- result will be a scalar
    else
      for i,v in ipairs(t) do
        if i ~= c then
          return v
        end
      end
    end

  end

  --
  -- Binary search
  --
  while (first <= last) do
    if tbl[middle][col] < q then
      first = middle + 1
    elseif tbl[middle][col] == q then
      break
    else
      last = middle - 1
    end
    middle = math.floor( (first + last) / 2 )
  end

  --
  -- no exact match
  --
  if first > last then
  
    -- linear interpolation
    if method == "linear" then
      local res = {}
      local fact = (q - tbl[last][col]) / (tbl[first][col] - tbl[last][col])
      local mincols = math.min(#tbl[first], #tbl[last])
      for i = 1, mincols do 
        if i ~= col then
          local val = tbl[last][i] + fact*(tbl[first][i]-tbl[last][i])
          table.insert(res, val)
        end
      end 
      return res
    -- nearest neighbor solution
    elseif method == "nearest" then
      local left = math.abs(q - tbl[last][col])
      local right = math.abs(tbl[first][col] - q)
      if left < right then
        return exclude_col(tbl[last], col)
      else
        return exclude_col(tbl[first], col)
      end
    -- return the next value within the bracket
    elseif method == "next" then
      return exclude_col(tbl[first], col)
    -- return the previous value within the bracket
    elseif method == "previous" then
      return exclude_col(tbl[last], col)
    -- error
    else
      error("Uknown interp method "..method)
    end

  --
  -- an exact match
  --
  else

    return exclude_col(tbl[middle], col)

  end

end
