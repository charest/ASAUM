-- http://lua-users.org/wiki/FileInputOutput
require "stringutil"

-- Find the length of a file
--   filename: file name
-- returns
--   len: length of file
--   asserts on error
function length_of_file(filename)
  local fh = assert(io.open(filename, "rb"))
  local len = assert(fh:seek("end"))
  fh:close()
  return len
end

-- Return true if file exists and is readable.
function file_exists(path)
  local file = io.open(path, "rb")
  if file then file:close() end
  return file ~= nil
end

-- Guarded seek.
-- Same as file:seek except throws string
-- on error.
-- Requires Lua 5.1.
function seek(fh, ...)
  assert(fh:seek(...))
end

-- Read an entire file.
-- Use "a" in Lua 5.3; "*a" in Lua 5.1 and 5.2
function readall(filename)
  local fh = assert(io.open(filename, "rb"))
  local contents = assert(fh:read(_VERSION <= "Lua 5.2" and "*a" or "a"))
  fh:close()
  return contents
end

-- Write a string to a file.
function write(filename, contents)
  local fh = assert(io.open(filename, "wb"))
  fh:write(contents)
  fh:flush()
  fh:close()
end

-- Read, process file contents, write.
function modify(filename, modify_func)
  local contents = readall(filename)
  contents = modify_func(contents)
  write(filename, contents)
end

-- Read csv files, takes 4 parameters.
-- path: the path of the CSV file to read - mandatory
-- sep: the separator character of the fields. Optionsl, defaults to ','
-- tonum: whether to convert fields to numbers if possible. Optional. Defaults to true
-- null: what value should null fields get. Optional. defaults to ''
function read(path, skiprows, sep, usecols, tonum, null)
  tonum = tonum or true
  sep = sep or ','
  null = null or ''
  skiprows = skiprows or 0

  if usecols and type(usecols) ~= "table" then
    usecols = {usecols}
  end

  local csvFile = {}
  local file = assert(io.open(path, "r"))
  local row = 0

  for line in file:lines() do
    if row >= skiprows then
      fields = line:split(sep)

      --
      -- Extract columnes
      --
      if usecols then
        -- single value
        if #usecols == 1 then
          fields = fields[usecols[1]]
        -- multiple columns
        else
          local newfields = {}
          for i,v in ipairs(usecols) do
            table.insert(newfields, fields[v])
          end
          fields = newfields
        end
      end

      --
      -- convert numeric fields to numbers
      --
      if tonum then
        -- numerious values
        if type(fields) == "table" then
          for i,v in ipairs(fields) do
            if v == '' then
              fields[i] = null
            end
            fields[i] = tonumber(v) or v
          end
        -- scalar value
        else
          fields = tonumber(fields) or fields
        end
      end

      table.insert(csvFile, fields)
    end
    row = row + 1
  end
  file:close()
  return csvFile
end

-- get the path of this script
function script_path()
  local str = debug.getinfo(2, "S").source:sub(2)
  return str:match("(.*/)") or "."
end

