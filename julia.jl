# activate project environment
using Pkg
Pkg.activate("C:\\Users\\Scuff\\OneDrive\\Course Documents\\Lobster Hatchery\\Code\\")

using RCall, DataFrames, Plots
plotlyjs()
theme(:ggplot2)
default(
    legend=:outertopright,
    legend_title="Legend",
    legend_background_color=:white,
    linewidth=3
    )

R"""
library(mizer)
library(mizerShelf)
"""
@rimport mizer

R"source('analysis/hatchery.R')"
function projectHatchery(params; t_max = 100, t_save=1, dt=0.1)
    rcall(R"projectHatchery", params, t_max = t_max, t_save = t_save, dt = dt)
end

function plotSpectra(object)
    biomass = mizer.getBiomass(object)
end

# define params we are using
params = R"NWMed_params";
params = mizer.steady(params);

biomass = plotSpectra(params);

# project params
sim = mizer.project(params);