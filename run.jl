using Cyton
using Gadfly
using CellCycle: createPopulation, stretchedCellFactory, BrduStimulus, dnaStainLevels, runModel!
struct RunResult
  plot::Union{Nothing, Plot}
  dnas::Vector{Real}
  brdus::Vector{Real}
  stimDur::Real
  model::CytonModel
  negDnaCnt::Int
end

forReal = true
results = RunResult[]
if forReal
  println("-------------------- start --------------------")
  stimDurs = [0.5, 1.0, 2.0, 4.0]
  for stimDur in stimDurs
    local model, rt, dnas, result, p, brdus, negDnaCnt
    model = createPopulation(200000, stretchedCellFactory)
    # model.properties[:Δt] = 0.01
  
    stim = BrduStimulus(0.5, 0.5+stimDur)
    rt = runModel!(model, 1.0+stimDur, stim)
  
    (dnas, brdus, negDnaCnt) = dnaStainLevels(model);

    title = "pulse=$(Int(round(stimDur*60)))mins"
    p = plot(x=dnas, 
    y=brdus, 
    Geom.histogram2d, 
    Guide.xlabel("DNA"), 
    Guide.ylabel("BrdU"), 
    Coord.cartesian(xmin=50000, xmax=200000, ymin=1, ymax=5), 
    Scale.y_log10, 
    Guide.title(title),
    Theme(background_color="white", key_position=:none))
    display(p)
    p |> PNG("/Users/thomas.e/Desktop/$(title).png", 15cm, 15cm)
    # p = nothing
    result = RunResult(p, dnas, brdus, stimDur, model, negDnaCnt)
    push!(results, result)
  end  
end

include("src/WriteFCS.jl")

for (i, result) in enumerate(results)
  dnas = result.dnas
  brdus = result.brdus
  stimDur = Int(round(result.stimDur*60))
  params = Dict(
    "\$TOT" => string(length(dnas)),
    "\$P1N" => "DNA",
    "\$P2N" => "BrdU",
    "\$PAR" => "2",
    "\$P1E" => "0,0",
    "\$P2E" => "0,0",
    "\$P1B" => string(sizeof(Float32)*8),
    "\$P2B" => string(sizeof(Float32)*8),
    "\$P1R" => "262144", #string(maximum(dnas)),
    "\$P2R" => "262144", #string(maximum(brdus)),
  )
  data = Dict(
    "BrdU" => convert(Vector{Float32}, brdus),
    "DNA" => convert(Vector{Float32}, dnas),
  )
  fcs = FlowSample{Float32}(data, params)
  writeFcs("/Users/thomas.e/Desktop/$(stimDur).fcs", fcs)
end

# tm = [r.stimDur*60 for r in results]
# negCnts = [r.negDnaCnt/length(r.model.agents)*100 for r in results]
# p = plot(x=tm, 
# y=negCnts,
# Guide.xlabel("Time (mins)"), 
# Guide.ylabel("%BrdU(-ve) DNA(x2)"), 
# Coord.cartesian(xmin=0, xmax=250, ymin=0, ymax=12), 
# Theme(background_color="white"))
# display(p)
# p |> PNG("/Users/thomas.e/Desktop/survival.png", 6inch, 6inch)
