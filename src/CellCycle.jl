module CellCycle

"""
Implements the Streched Cell Cycle Model:

Stretched cell cycle model for proliferating lymphocytes
Dowling, Mark R.Kan, AndreyHeinzel, SusanneZhou, Jie H. S.Marchingo, Julia M.
PROCEEDINGS OF THE NATIONAL ACADEMY OF SCIENCES 2014
https://dx.doi.org/10.1073/pnas.1322420111
"""

using Cyton
import Cyton: shouldDie, shouldDivide, inherit, step, stimulate
using Gadfly: plot, layer, cm, Gadfly, Theme, Guide, Geom, Col, mm, style, Scale, PNG, SVG, SVGJS, Coord, Plot
using Serialization, Cairo, DataFrames, Base.Threads, Colors, Cairo

# Gadfly defaults
Gadfly.set_default_plot_size(20cm, 20cm)
Gadfly.push_theme(Theme(background_color="white", alphas=[0.5]))

@enum Phase G1 S G2M
# CpG stimulated B cells
λ   = LogNormalParms(12.34, 3.48; natural=true)
kG1 = 0.27
kS  = 0.57
# BRDU labelled and unlabelled floursence level
brduLo = LogNormalParms(log(180), log(1.8))
brduHi = LogNormalParms(log(8000), log(1.8))
# 7AAD DNA 2x and 4x labelling
dnaLo = NormalParms(75000, 5000)
dnaHi = NormalParms(150000, 10000)

mutable struct BrduStatus
  lo::Float64
  hi::Float64
  current::Float64
  positive::Bool
end
function BrduStatus(brduLo::DistributionParmSet, brduHi::DistributionParmSet)
  l = Cyton.draw(brduLo)
  h = Cyton.draw(brduHi)
  BrduStatus(l, h, l, false)
end

mutable struct CycleTimer <: FateTimer
  divTime::Time
  kG1::Float64
  kS::Float64
  startTime::Time
  brdu::BrduStatus
end

"Constructor for new cells at t=0"
function CycleTimer(λ::DistributionParmSet, 
  kG1::Float64, 
  kS::Float64,
  brduLo::DistributionParmSet, 
  brduHi::DistributionParmSet)
  
  # Timers are desynchronised by distributing them randomly 
  # over their cycle.
  l = Cyton.draw(λ)
  s = - l * rand()
  
  CycleTimer(l, kG1, kS, s, BrduStatus(brduLo, brduHi))
end

function phase(cycle::CycleTimer, time::Time)
  "Returns the phase of the cycle and the proportion of time spent in that phase"
  δt = time - cycle.startTime
  divTime = cycle.divTime
  
  kG1 = cycle.kG1 * divTime
  kS = cycle.kS * divTime
  if δt <= kG1
    return (G1, δt/kG1)
  end
  if δt <= kG1 + kS
    return (S, (δt-kG1)/kS)
  end
  return (G2M, (δt-kG1-kS)/(kG1+kS))
end

"The step function for the cycle timer"
function step(cycle::CycleTimer, time::Time, Δt::Duration)
  if time >= cycle.divTime + cycle.startTime
    # Cell has (fake) divided, reset the timer and BrdU level.
    cycle.startTime = time
    # cycle.brdu.current = cycle.brdu.lo
    # cycle.brdu.positive = false
  end
  return nothing
end

remaining(cycle::CycleTimer, time::Time) = cycle.divTime + cycle.startTime - time
remaining(cell::Cell, time::Time) = remaining(cell.timers[1], time)

function stretchedCellFactory(birth::Float64=0.0)
  cell = Cell(birth)
  addTimer(cell, CycleTimer(λ, kG1, kS, brduLo, brduHi))
  return cell
end

function runModel!(model::CellPopulation, runDuration::Float64, stimulus::Stimulus, callback::Function=(_) -> nothing)
  Δt = modelTimeStep(model)
  for _ in 0:Δt:runDuration
    step(model, stimulus)
    callback(model)
  end
end

function survivalCurves(model::CellPopulation)
  totalAlpha = Vector{Float64}()
  g1Alpha = Vector{Float64}()
  sg2mAlpha = Vector{Float64}()
  
  runDuration = model_time(model)
  for cell in model.cells
    push!(totalAlpha, remaining(cell, runDuration))
    cycle = cell.timers[1]
    if phase(cycle, runDuration) == G1
      r = cycle.divTime * cycle.kG1 + cycle.startTime - runDuration
      push!(g1Alpha, r)
    else
      r = (cycle.divTime + cycle.startTime - runDuration) * (1 - cycle.kG1)
      push!(sg2mAlpha, r)
    end
  end

  sort!(totalAlpha)
  n = length(totalAlpha)
  tmT = 1 .- collect(1:1:n) ./ n
  sort!(g1Alpha)
  n = length(g1Alpha)
  tmG1 = 1 .- collect(1:n)/n
  sort!(sg2mAlpha)
  n = length(sg2mAlpha)
  tmSg2m = 1 .- collect(1:n)/n

  h = plot()
  plot!(g1Alpha, tmG1, label="G1", lw=3, lc="green")
  plot!(sg2mAlpha, tmSg2m, label="S/G2/M", lw=3, lc="red")
  plot!(totalAlpha, tmT, label="Total", lw=3, lc="blue")
  xlabel!("Time (h)")
  ylabel!("proportion")
  yaxis!(:log)
  yticks!([10.0^x for x in [-3, -2, -1, 0]])
  display(h)

  println("Done at model time: $(model.properties[:step_cnt]*model.properties[:Δt])")
end

function dnaStainLevels(model::CellPopulation)
  time = modelTime(model)
  NCells = length(model.cells)
  brdus = zeros(NCells)
  dnas = zeros(NCells)
  negDnaCnt = 0
  for (i, c) in enumerate(model.cells)
    fate = c.timers[1]
    brdus[i] = fate.brdu.current
    (p, timeInPhase) = phase(fate, time)
    
    l = Cyton.draw(dnaLo)
    if p == G1
      dnas[i] = l
      continue
    end

    h = Cyton.draw(dnaHi)
    if p == S
      t = timeInPhase
      dnas[i] = l + (h-l)*t
      continue
    end
    
    if p == G2M
      dnas[i] = h
      if !fate.brdu.positive
        negDnaCnt += 1
      end
    end
  end
  return (dnas, brdus, negDnaCnt)
end

struct BrduStimulus <: Stimulus
  pulseStart::Real
  pulseEnd::Real
end

function stimulate(cell::Cell, stim::BrduStimulus, time::Time)

  pulseStart = stim.pulseStart
  pulseEnd   = stim.pulseEnd
  inPulse = pulseStart ≤ time ≤ pulseEnd
  if !inPulse
    return
  end

  cycle = cell.timers[1]
  (p, timeInS) = phase(cycle, time)
  if p ≠ S
    return
  end

  brdu = cycle.brdu
  brdu.positive = true
  brdu.current = brdu.lo + timeInS * (brdu.hi - brdu.lo)

end

end # module CellCycle
