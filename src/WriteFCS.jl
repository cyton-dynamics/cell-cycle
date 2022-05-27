# module WriteFCS
# export writeFcs

using FileIO, CRC32c, FCSFiles

const SEP = 0x0C
const HEADER_LENGTH = 256

function get_endian_func(byte_order::String)
  # Taken from FCSFiles repo
  if byte_order == "1,2,3,4" # least significant byte first
      return htol
  elseif byte_order == "4,3,2,1" # most significant byte first
      return hton
  else
      error("Byte order is not supported.")
  end
end

struct FcsOffsets
  textStart::Int
  textEnd::Int
  dataStart::Int
  dataEnd::Int
  analysisStart::Int
  analysisEnd::Int
  function FcsOffsets(;textStart=0, textEnd=0, dataStart=0, dataEnd=0, analysisStart=0, analysisEnd=0)
    new(textStart, textEnd, dataStart, dataEnd, analysisStart, analysisEnd)
  end
end

function fcsOffsets(params::Dict{String, String}, data::Dict{String, Vector{T}}) where T <: AbstractFloat
  textStart = 256
  textEnd = textStart
  for (k, v) in params
    textEnd += 2 + length(k) + length(v)
  end
  dataStart = textEnd + 1
  dataEnd = dataStart 
  for (_, v) in data
    dataEnd += length(v) * sizeof(T)
  end

  FcsOffsets(textStart=textStart, textEnd=textEnd, dataStart=dataStart, dataEnd=dataEnd)
end

function fcsHeader(offsets::FcsOffsets)
  h = rpad("FCS3.0", 9, " ") # Why 9?
  h *= lpad(string(offsets.textStart), 8, " ")
  h *= lpad(string(offsets.textEnd), 8, " ")
  h *= lpad(string(offsets.dataStart), 8, " ")
  h *= lpad(string(offsets.dataEnd), 8, " ")
  h *= lpad(string(offsets.analysisStart), 8, " ")
  h *= lpad(string(offsets.analysisEnd), 8, " ")
  return rpad(h, HEADER_LENGTH, " ")
end

function fcsText(params::Dict{String, String}) where T <: AbstractFloat
  sb = IOBuffer()

  for (k, v) in params
    write(sb, SEP)
    write(sb, k)
    write(sb, SEP)
    write(sb, v)
  end
  write(sb, SEP)

  return take!(sb)
end

function fcsData(params::Dict{String, String}, data::Dict{String, Vector{T}}) where T <: AbstractFloat
  eventCnt = parse(Int64, params["\$TOT"])
  paramCnt = parse(Int64, params["\$PAR"])
  sb = IOBuffer()

  byte_mapping = get_endian_func(params["\$BYTEORD"])

  for e in 1:eventCnt
    for p in 1:paramCnt
      parm = "\$P" * string(p) * "N"
      key = params[parm]
      d = data[key][e]
      b = byte_mapping(d)
      write(sb, b)
    end
  end

  return take!(sb)
end

function writeFcs(fcsFile::String, fcs::FCSFiles.FlowSample{T}) where T <: AbstractFloat
  fh = open(fcsFile, "w")
  try
    writeFcs(fh, fcs)
  finally
    close(fh)
  end
end

function fixParams(params::Dict{String, String})
  blanks = lpad("", 8, " ")
  required_params = Dict{String, String}(
    "\$BEGINANALYSIS" => "0", # Byte-offset to the beginning of the ANALYSIS segment.
    "\$BYTEORD" => "4,3,2,1", # Byte order for data acquisition computer.
    "\$DATATYPE" => "F", # Type of data in DATA segment (ASCII, integer, floating point).
    "\$ENDANALYSIS" => "0", # Byte-offset to the last byte of the ANALYSIS segment.
    "\$MODE" => "L", # Data mode (list mode - preferred, histogram - deprecated).
    "\$NEXTDATA" => "0", # Byte offset to next data set in the file.
    "\$BEGINDATA" => blanks, # Byte-offset to the beginning of the DATA segment.
    "\$BEGINSTEXT" => "0", # Byte-offset to the beginning of a supplemental TEXT segment.
    "\$ENDDATA" => blanks, # Byte-offset to the last byte of the DATA segment.
    "\$ENDSTEXT" => "0", # Byte-offset to the last byte of a supplemental TEXT segment.
    )

  fixed = copy(params)
  for (k, v) in required_params
    if !haskey(fixed, k)
      fixed[k] = v
    end
  end

  return fixed
end

function writeFcs(fcsFile::IOStream, fcs::FCSFiles.FlowSample{T}) where T <: AbstractFloat
  params = fixParams(fcs.params)
  data = fcs.data

  offsets = fcsOffsets(params, data)
  println(offsets)

  params["\$BEGINDATA"] = lpad(string(offsets.dataStart), 8, " ")
  params["\$ENDDATA"] = lpad(string(offsets.dataEnd), 8, " ")

  crc::UInt32 = 0x00000000
  l = 0

  txt = fcsHeader(offsets)
  crc = crc32c(txt, crc)
  l += write(fcsFile, txt)
  
  txt = fcsText(params)
  crc = crc32c(txt, crc)
  l += write(fcsFile, txt)
  
  txt = fcsData(params, data)
  crc = crc32c(txt, crc)
  l += write(fcsFile, txt)
    
  l += write(fcsFile, crc)

  return l
end

