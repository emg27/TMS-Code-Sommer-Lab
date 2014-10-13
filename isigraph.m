function[N,isi]=isigraph(t_clust,t_min,t_max,bin,last_bin)
%%ISI Values
isi_range=find(t_clust>t_min & t_clust<=t_max);

isi=diff(t_clust(isi_range));
if size(isi)>0
edges=0:bin:last_bin;
[N,isi_bin]=histc(isi,edges);
bar(edges,N,'histc');
end
