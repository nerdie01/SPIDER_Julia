module VizUtils

using GLMakie

function quickmap(F::Matrix)
    fig, ax, hm = heatmap(real.(F))
    Colorbar(fig[:, end+1], hm)

    fig
end

function quickanim(F::Array, steps::Int, fps::Int=20, name::String="quickanim.gif")
    fig, _, hm = heatmap(real.(F[:,:,1]))
    Colorbar(fig[:, end+1], hm)
    record(fig, name, framerate=fps, 1:steps) do i
        hm[1] = F[:,:,i]
    end
end

end