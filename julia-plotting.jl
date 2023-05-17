# activate project environment
using Pkg
Pkg.activate("C:\\Users\\Scuff\\OneDrive\\Course Documents\\Lobster Hatchery\\Code\\")

# import packages
using RCall, DataFrames, Plots

# color palettes
ibm = palette(["#6a90ff", "#7660f2", "#d22680", "#f46000", "#f8af00"]) # ibm colorblind palette
wong = palette(["#e09e00", "#65b5ea", "#359e72", "#efe337", "#2473b3", "#cd5d00", "#c579a8", "#000000"])
tol = palette(["#322489", "#2b7731", "#56aa99", "#91ccef", "#dbcb74", "#c46677", "#a3459a", "#822256"])

# plot theme/default styling
theme(:dao) # nicer style
# theme(:ggplot2) # for consistency with mizer Plots
# plotlyjs() # backend for consistency with mizer plots
default(
    # legend
    legend=:outertopright,
    legend_title="Legend",
    legend_background_color=:white, # for ggplot2 theme

    # sizing
    size = [800,360], # aspect ratio
    dpi=300, # high quality
    margin=0.6Plots.cm, # keep text in plot frame
    widen=true, # automatically widen axis limits
    
    # colors/line styling
    palette = wong, # colorblind_palette
    linewidth=2.4,
    
    # fonts
    fontfamily="Computer Modern", # LaTeX font
    guidefontsize = 14,
    tickfontsize = 12,
    legendfontsize = 10,
    # legendtitlefont = 12
    )

# Rcall imports
R"""
library(mizer)
library(mizerShelf)
"""

@rimport mizer
@rimport base as R

R"source('analysis/hatchery.R')"

# Functions to interface with mizer + plotting
function projectHatchery(params; t_max = 100, t_save=1, dt=0.1)
    rcall(R"projectHatchery", params, t_max = t_max, t_save = t_save, dt = dt)
end

issim(object) = rcopy(R.class(object)[1]) == "MizerSim"
function getSpecies(object)
    if issim(object)
        params = object["params"]
    else
        params = object
    end

    return rcopy(R.rownames(mizer.species_params(params)))
end

function getN(object)
    if issim(object)
        N = rcopy(object["n"])
        t = parse.(Float64, rcopy(R.rownames(object["n"])))
    else
        N = rcopy(object["initial_n"])
        t = [0.0]
    end

    return t, N
end
function getBiomass(object)
    rN = mizer.getBiomass(object)
    rt = R.rownames(rN)
    
    t = rcopy(rt)
    if isnothing(t)
        t = [0.0]
    else
        t = parse.(Float64, t)
    end

    N = rcopy(rN)

    return t, N
end
function getYield(object)
    rN = mizer.getYield(object)
    rt = R.rownames(rN)
    
    t = rcopy(rt)
    if isnothing(t)
        t = [0.0]
    else
        t = parse.(Float64, t)
    end

    N = rcopy(rN)

    return t, N
end

function plotBiomassRelative(
    sim;
    species = getSpecies(sim), tspan=nothing, kwargs...)
    
    @assert issim(sim)
    @assert issubset(species, getSpecies(sim))
    
    t, N = getBiomass(sim)

    initial_N = N[1, :]
    for i = axes(N, 1)
        N[i, :] = 100 * (N[i, :] - initial_N) ./ initial_N
    end

    if isnothing(tspan)
        t_indices = ones(Bool, length(t))
    else
        t_indices = map(t -> tspan[1] ≤ t ≤ tspan[2], t)
    end

    plot(
        xaxis = "Time [years]",
        yaxis = "Biomass change %",
        kwargs...
        )
    for (i, s) in enumerate(getSpecies(sim))
        if s in species
            plot!(t[t_indices], N[t_indices, i], label = s)
        end
    end
    plot!()
end

function plotYieldRelative(
    sim;
    species = getSpecies(sim), tspan=nothing, kwargs...)
    
    @assert issim(sim)
    @assert issubset(species, getSpecies(sim))
    
    t, N = getYield(sim)

    initial_N = N[1, :]
    for i = axes(N, 1)
        N[i, :] = 100 * (N[i, :] - initial_N) ./ initial_N
    end

    if isnothing(tspan)
        t_indices = ones(Bool, length(t))
    else
        t_indices = map(t -> tspan[1] ≤ t ≤ tspan[2], t)
    end


    plot(
        xaxis = "Time [years]",
        yaxis = "Yield change %",
        kwargs...
        )
    for (i, s) in enumerate(getSpecies(sim))
        if s in species
            plot!(t[t_indices], N[t_indices, i], label = s)
        end
    end
    plot!()
end

function plotSpectra(object;
    species = getSpecies(object), t = nothing,
    kwargs...)

    species = filter(s -> s in species, getSpecies(object))
    species_indices = map(s -> s in species, getSpecies(sim))

    t_arr, N = getN(object)
    w = rcopy(issim(object) ? object["params"]["w"] : object["w"])
    
    if isnothing(t)
        t = maximum(t_arr)
    else
        @assert t ∈ t_arr
    end

    t_index = findfirst(x -> x == t, t_arr)
    
    # filter N
    N = issim(object) ? N[t_index, species_indices, :] : N[species_indices, :]

    # compute biomass
    biomass = zeros(size(N))
    for i in axes(biomass, 1)
        biomass[i,:] .= N[i,:] .* w
    end
    y_min = 1e-7
    biomass[biomass .<= y_min] .= 0.0

    show(biomass)
    notzero(x) = !iszero(x)

    plot(
        scale = :log,
        xaxis = "Size [g]",
        yaxis = "Biomass density";
        kwargs...)
    for i in axes(biomass,1)
        w_indices = notzero.(biomass[i, :])
        plot!(w[w_indices], biomass[i, w_indices], label = species[i])
    end

    plot!()
end

# example usage ----------------------------------------------------------------
params = R"NWMed_params";
params = mizer.steady(params);

# project params
sim = mizer.project(params, t_save = 0.2, t_max = 30, effort = 10);

# plotting 
target = ["Angler fish", "Blue whiting",
"Horned octopus", "Shortfin squid", "Angular crab", "Harbour crab"]

# plot size spectrum
plotSpectra(sim, species = target, ylims = (1e-7, :auto), t = 0.0)
plotSpectra(sim, species = target, ylims = (1e-7, :auto))

# plot biomass
plotBiomassRelative(sim, species = target, tspan = (0, 20))

# plot yield
plotYieldRelative(sim, species = target, tspan = (0,20))

# save last figure created
# savefig("analysis/test.png")

# multispecies experiment conditions: ------------------------------------------

# determine more reasonable hatchery input by fraction of total number of mature
# individuals maturity
N_mat = mizer.getN(params, min_w=mizer.species_params(params)[:w_mat])

# number of individuals to add from hatchery as a % of N_mat, i.e. ratio=1 => 
# add one individual for every mature individual in population
ratio = 1

hatchery_target = "Angular crab"
hatchery_params = DataFrame(species = hatchery_target, mu = 0.5, sigma = 0.1, annual_N=N_mat[hatchery_target]*ratio)
@rput hatchery_params

# define new hatchery time dependence function (default is constant 1)
# R"""
# hatchery_T <- function(t) {
#     # return(1 - sign(sin(2 * pi * t))) # square wave
#     return(1)
# }
# """;

baseline = mizer.project(params, t_max=30, t_save=0.2, dt=0.05);
hatchery = projectHatchery(params, t_max=30, t_save=0.2, dt=0.05);

plot_baseline_spectra = plotSpectra(baseline, species = target, ylims = (1e-7, :auto))
plot_hatchery_spectra = plotSpectra(hatchery, species = target, ylims = (1e-7, :auto))

plot_biomass = plotBiomassRelative(hatchery, species=target, tspan = (0, 15))
plot_yield = plotYieldRelative(hatchery, species=target, tspan = (0, 15))

display.([plot_baseline_spectra, plot_hatchery_spectra, plot_biomass, plot_yield])

# Comparable fishing effort calculations ---------------------------------------

function shannon(N)
    n = filter(!iszero, N) # remove zeros (upsets log function)

    total = sum(n)
    p = n / total
    
    return - sum(p .* log.(p))
end

t = getN(baseline)[1]

N_baseline = rcopy(mizer.getN(baseline))
H_baseline = shannon.(eachrow(N_baseline))

N_hatchery = rcopy(mizer.getN(hatchery))
H_hatchery = shannon.(eachrow(N_hatchery))


H_baseline = shannon.(eachrow(N_baseline))
H_hatchery = shannon.(eachrow(N_hatchery))

H_rel = (H_hatchery .- H_baseline) ./ H_baseline * 100

plot(
    t,
    H_rel,
    label="Hatchery",
    yaxis = "Relative Shannon index %",
    xaxis = "Years",
    legend = :none,
    size = [600,360],
    # margin = 0.5Plots.cm,
    dpi=300,
    linewidth=3
)
savefig("Shannon3.png")

function shannonfromsim(params; kwargs...)

    sim = mizer.project(params; kwargs...)

    N_sim = rcopy(mizer.getN(sim))
    H_sim = shannon(N_sim[end,:]) # Shannon index for final row

    return (H_sim - H_baseline[end]) / H_baseline[end] * 100
end

for e ∈ [1,2,3,4,5] # 1 is the `initial_effort` (default) for these params
    @info e, shannonfromsim(params; t_max = 30, t_save=0.2, dt=0.05, effort = e)
end
# suggests that hatchery with 1x mature N is as bad for ecosystem as ~4x fishing
# effort - probably not sensible to use since fishing initially increases
# Shannon index

begin
    H_avg = 0.0
    for i = 1:25
        local N = copy(N_baseline[end, :])
        N[i] /= 1.369
        H_avg += shannon(N)
    end
    H_avg /= 25
    (H_avg - H_baseline[end]) / H_baseline[end] * 100
end

# average change to H from extinction of a species in the model is -3.27% to put
# my numbers into context

# average change to H from halving the population of a single species in the
# model is -0.357%

H_rel[end]

# species checks ---------------------------------------------------------------

params2 = params;
rcopy(mizer.species_params(params2))
new_n = rcopy(params2["initial_n"])
new_n[25,:] .= 0.0
params2["initial_n"] = new_n

params2 = mizer.steady(params2);
sim = projectHatchery(params2, t_max=30, t_save=0.2, dt=0.05);

plotBiomassRelative(sim, species=target, tspan = (0, 15))