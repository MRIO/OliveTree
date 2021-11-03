function cmap1 = gen_divergent_colormap()
% colormap generated 
s=0:1/256:1-(1/256);
rgb1 =[0 0 1];
rgb2=[1 0 0];
cmap1 =  diverging_map(s,rgb1,rgb2);
colormap(cmap1);