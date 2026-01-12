module VizUtils

using GLMakie

function quickmap(F::Matrix, steps)
    fig, ax, hm = heatmap(real.(F))
    Colorbar(fig[:, end+1], hm)

    fig
end

function quickanim(F::Array, steps)
    fig, _, hm = heatmap(real.(F[:,:,1]))
    Colorbar(fig[:, end+1], hm)
    record(fig, "quickanim.mp4", 1:steps) do i
        hm[1] = F[:,:,i]
    end
end

end