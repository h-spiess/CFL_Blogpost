#%%
using Distributions
using Random
using CairoMakie, AlgebraOfGraphics
set_aog_theme!()

#%%
function sample_eda(lat, hue)
    red = @. 0 < hue <= 90 || 270 < hue <= 360

    eda = map(zip(lat, red)) do (l, r)
        if r && (l <= 45)
            rand(Uniform(0.5, 1.))
        elseif !r && (l <= 45)
            b = rand(Bernoulli(2/3))
            b*rand(Uniform(0.0, 0.5)) + (1-b)*rand(Uniform(0.5, 1.0))
        elseif !r && (l > 45)
            rand(Uniform(0.0, 0.5))
        elseif r && (l > 45)
            b = rand(Bernoulli(1/5))
            b*rand(Uniform(0.0, 0.5)) + (1-b)*rand(Uniform(0.5, 1.0))
        else
            nothing
        end
    end
    eda
end

function sample_hue(lat)
    hue = map(lat) do l
        if l <= 45
            rand(Uniform(0, 180))
        else 
            rand(Uniform(180, 360))
        end
    end
    hue
end

function sample_hue_eda(;nsamples=2000)
    # hue -> eda (electrodermal response) / skin conductance

    # latent unobserved confounder lat(itude)
    lat = rand(Uniform(0, 90), nsamples)

    hue = sample_hue(lat)

    eda = sample_eda(lat, hue)
    (hue, eda, lat)
end

function sample_eda_from_manipulation(man_hue; nsamples=100)
    # man(hue) -> eda (electrodermal response) / skin conductance

    # latent unobserved confounder lat(itude)
    lat = rand(Uniform(0, 90), nsamples)

    hue = repeat([man_hue], nsamples)

    eda = sample_eda(lat, hue)
    (hue, eda, lat)
end

function sample_hue_eda_unconfounded(; nsamples=2000)
    # hue -> eda (electrodermal response) / skin conductance

    # latent unobserved confounder lat(itude)
    lat = rand(Uniform(0, 90), nsamples)

    hue = rand(Uniform(0, 360), nsamples)

    red = @. 0 < hue <= 90 || 270 < hue <= 360

    eda = sample_eda(lat, red)
    (hue, eda, lat)
end

function densityfn(h, e)
    @assert 0 <= h <= 360 "Invalid value for hue: $h"
    @assert 0 <= e <= 1.0 "Invalid value for eda: $e"

    if (0 <= h <= 90) && (e >= 0.5)
        1.0*2
    elseif (90 < h <= 180) && (e < 0.5)
        2/3*2
    elseif (90 < h <= 180) && (e >= 0.5)
        1/3*2
    elseif (180 < h <= 270) && (e < 0.5)
        1.0*2
    elseif (270 < h <= 360) && (e < 0.5)
        1/5*2
    elseif (270 < h <= 360) && (e >= 0.5)
        4/5*2
    else
        0.0
    end
end

#%%
function sample_density_plot()
    hue, eda, _ = sample_hue_eda()

    sample_density_plot(hue, eda; title="Samples and true density")
end

function sample_density_plot(hue, eda; edaname="eda", densityfunc=densityfn, densityname="Cond. density p(eda|hue)", title="", dotcolor=:black, markersize=9)
    h = range(0,360,length=200)  # note ': this is a row vector
    e = range(0,1.,length=200)'
    d = @. float(densityfunc(h, e))

    df = (hue=hue, eda=eda)
    hueeda = mapping([h] => "hue", [e'] => edaname, [d] => densityname) * visual(Heatmap, colormap=:Blues_5) + data(df) * mapping(:hue, :eda => edaname) * visual(Scatter, color=dotcolor, markersize=markersize)
    draw(hueeda, axis=(
        xticks=[0., 90., 180., 270., 360.],
        yticks=[0., 0.5, 1.0],
        title=title
    ))
end

#%%
function estimate_density_by_counts(hue, eda)
    buckets = zeros(4, 2)
    for (h, e) in zip(hue, eda)
        if (0 <= h <= 90) && (e < 0.5)
            buckets[1, 1] += 1
        elseif (0 <= h <= 90) && (e >= 0.5)
            buckets[1, 2] += 1
        elseif (90 < h <= 180) && (e < 0.5)
            buckets[2, 1] += 1
        elseif (90 < h <= 180) && (e >= 0.5)
            buckets[2, 2] += 1
        elseif (180 < h <= 270) && (e < 0.5)
            buckets[3, 1] += 1
        elseif (180 < h <= 270) && (e >= 0.5)
            buckets[3, 2] += 1
        elseif (270 < h <= 360) && (e < 0.5)
            buckets[4, 1] += 1
        elseif (270 < h <= 360) && (e >= 0.5)
            buckets[4, 2] += 1
        else
            error("Invalid hue and eda: $h, $e.")
        end
    end
    buckets ./ sum(buckets; dims=2)
end

#%%
function estimate_density_by_counts_by_color(hue, eda)
    buckets = zeros(2, 2)
    red = @. 0 < hue <= 90 || 270 < hue <= 360
    
    for (r, e) in zip(red, eda)
        if r && (e < 0.5)
            buckets[1, 1] += 1
        elseif r && (e >= 0.5)
            buckets[1, 2] += 1
        elseif !r && (e < 0.5)
            buckets[2, 1] += 1
        elseif !r && (e >= 0.5)
            buckets[2, 2] += 1
        else
            error("Invalid hue and eda: $h, $e.")
        end
    end
    buckets ./ sum(buckets; dims=2)
end

# #%%
# hue, eda, _ = sample_hue_eda(;nsamples=200000)
# estimate_density_by_counts(hue, eda)

# hue, eda, _ = sample_eda_from_manipulation(200.; nsamples=200000)
# estimate_density_by_counts(hue, eda)

# #%%
# hue, eda, _ = sample_hue_eda_unconfounded(;nsamples=1000)
# sample_density_plot(hue, eda)

# hue, eda, _ = sample_hue_eda_unconfounded(;nsamples=20000)
# estimate_density_by_counts(hue, eda)
# estimate_density_by_counts_by_color(hue, eda)