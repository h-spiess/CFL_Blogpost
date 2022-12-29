using AlgebraOfGraphics, CairoMakie
CairoMakie.activate!()
set_aog_theme!()
using Distributions
using Random
Random.seed!(42) #for reproducibility
using StatsBase
using Clustering
using KernelDensity
using Distances

include("hue_eda_data.jl")

hue, eda, _ = sample_hue_eda()

# reshaping into a form, which is more convenient to work with
hue, eda = reshape(hue, 1, :), reshape(eda, 1, :)

sample_density_plot(hue[1, :], eda[1, :]; title="Samples and true density")


#joint and marginal probability density estimators
kde_hueeda = kde([hue; eda]')
kde_hue = kde(hue[1, :]);


kde_hueeda = InterpKDE(kde_hueeda) #makes the prediction faster
kde_hue = InterpKDE(kde_hue);   #makes the prediction faster


pdf_eda_given_hue(hue, eda) = pdf(kde_hueeda, hue, eda)/pdf(kde_hue, hue);


sample_density_plot(hue[1, :], eda[1, :]; densityfunc=pdf_eda_given_hue,
    densityname="Estimated cond. density p(eda|hue)", 
    title="Samples and estimated density")


edarange = 0.0:0.0005:1.0
huerange = 0.0:0.025:360.0

gridpoints = [(h, e) for e in edarange, h in huerange]

grid_model_dens = [pdf_eda_given_hue(p...)
                   for p in gridpoints]

# coarsening for visualization
coarser_gridpoints = gridpoints[begin:20:end, begin:500:end]

sample_density_plot(first.(reshape(coarser_gridpoints, :)), last.(reshape(coarser_gridpoints, :)); 
    densityfunc=pdf_eda_given_hue,
    densityname="Estimated cond. density p(eda|hue)", 
    title="Grid points on the estimated density", dotcolor=:purple, markersize=2.5)


kmhue = kmeans(grid_model_dens, 4);


function kmpredict(kmresult, X)
    dmat = Distances.pairwise(Distances.Euclidean(), kmresult.centers, X)
    getindex.(argmin(dmat, dims=1), 1)[1, :]
end;


grid_model_dens_for_hue(hue) = [pdf_eda_given_hue(h, e)
                                for e in edarange, h in hue];


function sample_cluster_plot(hue, eda, assigned_clusters; title="")
    df = (hue=hue, eda=eda, Cluster=assigned_clusters)
    hueeda = data(df) * mapping(:hue, :eda, color=:Cluster) * visual(Scatter)

    nclusters = length(unique(assigned_clusters))
    colors = cgrad(:batlow, nclusters, categorical=true)
    draw(hueeda, axis=(
        xticks=[0., 90., 180., 270., 360.],
        yticks=[0., 0.5, 1.0],
        title=title,
    ), palettes=(color=colors,))
end;


assigned_clusters = Symbol.(kmpredict(kmhue,
    grid_model_dens_for_hue(hue[1, :])))
sample_cluster_plot(hue[1, :], eda[1, :], assigned_clusters;
    title="Samples clustered by hue\n(Observational partition of hue)")


kmeda = kmeans(grid_model_dens', 2)

# coarsening for visualization
coarser_gridpoints = gridpoints[begin:100:end, begin:100:end]

sample_density_plot(first.(reshape(coarser_gridpoints, :)), last.(reshape(coarser_gridpoints, :)); 
    densityfunc=pdf_eda_given_hue,
    densityname="Estimated cond. density p(eda|hue)", 
    title="Grid points on the estimated density", dotcolor=:purple, markersize=2.5)


grid_model_dens_for_eda(eda) = [pdf_eda_given_hue(h, e)
                                for e in eda, h in huerange]';


assigned_clusters = Symbol.(kmpredict(kmeda,
    grid_model_dens_for_eda(eda[1, :])))
sample_cluster_plot(hue[1, :], eda[1, :], assigned_clusters; 
    title="Samples clustered by eda\n(Observational partition of eda)")


hue_experiments = [50., 140., 220., 340.];


man_joint = map(hue_experiments) do hue_exp
    man_hue, eda_exp, _ = sample_eda_from_manipulation(hue_exp;
        nsamples=30)
    (man_hue, eda_exp)
end

man_hues = reduce(vcat, first.(man_joint))
eda_exps = reduce(vcat, last.(man_joint));


function manipulation_plot(man_hues, eda_exps)
    df = (hue=man_hues, eda=eda_exps)
    hueeda = data(df) * mapping(:hue, :eda) * visual(Scatter)
    draw(hueeda; axis=(
        xticks=[0., 90., 180., 270., 360.],
        yticks=[0., 0.5, 1.0],
        title="Experimental eda responses to manipulations"
    ))
end

manipulation_plot(man_hues, eda_exps)


assigned_E = kmpredict(kmeda, grid_model_dens_for_eda(eda_exps))
assigned_H = kmpredict(kmhue, grid_model_dens_for_hue(man_hues));


function conditional_probabilities_from_counts(assigned_H, assigned_E)
    buckets = map(Iterators.product(sort(unique(assigned_H)),
        sort(unique(assigned_E)))) do (h, e)
            sum(assigned_E[assigned_H.==h] .== e)
    end
    buckets ./ sum(buckets; dims=2)
end

p_conds = conditional_probabilities_from_counts(assigned_H, assigned_E)


kmjoint = kmeans([hue; eda], 8)

assigned_clusters = Symbol.(assignments(kmjoint))
sample_cluster_plot(hue[1, :], eda[1, :], assigned_clusters;
    title="Joint clustering")


dHue = fit(ZScoreTransform, hue)
hue_standardized = StatsBase.transform(dHue, hue)

dEda = fit(ZScoreTransform, eda)
eda_standardized = StatsBase.transform(dEda, eda)

kmjoint = kmeans([hue_standardized; eda_standardized], 8)

assigned_clusters = Symbol.(assignments(kmjoint))
sample_cluster_plot(hue[1, :], eda[1, :], assigned_clusters;
    title="Joint clustering (standardized)")


kmjoint = kmeans([hue_standardized; eda_standardized], 7)

assigned_clusters = Symbol.(assignments(kmjoint))
sample_cluster_plot(hue[1, :], eda[1, :], assigned_clusters;
    title="Joint clustering (standardized)")

