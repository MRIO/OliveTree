function [f2 f1] = show_cells(info)

[xx yy] = get_coords(info);

    function plot_via_scatter(vec)
        scatter(xx,yy,50,vec,'s','filled')
        axis tight
        axis off
    end

    function plot_via_pcolor(vec)
        vec = vec(:);
        xs = xx-min(xx)+1;
        ys = yy-min(yy)+1;
        mx = max(xs);
        my = max(ys);
        p = zeros(mx,my)+nan;
        for i=1:numel(xx)
            p(xs(i), ys(i)) = vec(i);
        end
        pcolor(p'), shading flat, axis tight, axis off
    end

f2 = @plot_via_scatter;
f1 = @plot_via_pcolor;
end