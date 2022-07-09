-- split a string
function string:split(sep, plain, max)
  sep = sep or ','

  if max and max < 1 then
    max = nil
  end

  local aRecord = {}

  if self:len() > 0 then
    max = max or -1

    local nField, nStart = 1, 1
    local nFirst,nLast = self:find(sep, nStart, plain)
    
    while nFirst and max ~= 0 do
      aRecord[nField] = self:sub(nStart, nFirst-1)
      nField = nField+1
      nStart = nLast+1
      nFirst,nLast = self:find(sep, nStart, plain)
      max = max-1
    end
    aRecord[nField] = self:sub(nStart)
  end

  return aRecord
end 
