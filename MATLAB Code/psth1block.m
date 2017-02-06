function[spk_d,trl_fr,bin_start_times,baseline,mean_trl_fr,binned_spks]=psth1block(pulse,time_b,time_a,cluster,gauss_size,b_correct)

bin_size=1;%10;
gauss_size=gauss_size;%/1000;
gau_time_vec=-2*gauss_size:1:2*gauss_size;% -2*gauss_size:.001:2*gauss_size;
sigma=gauss_size; mu=0;
gaus_ker1 = normpdf(gau_time_vec,mu,sigma).*bin_size;
gaus_ker=gaus_ker1/sum(gaus_ker1);

trl_fr=[];
for n=1:length(pulse)
    TMS_art=pulse(n);
    place=find(cluster>(TMS_art-time_b) & cluster<(TMS_art+time_a));
    spks=cluster(place)-TMS_art;

   % bin_size=1;
    bin_start_times=-time_b:bin_size:time_a;
    binned_spks=zeros(length(bin_start_times),1);
    
    for bin_i=1:length(bin_start_times)
        binned_spks(bin_i)=sum(spks>=(-time_b+bin_size*(bin_i-1)) &...
            spks<(-time_b+bin_size*bin_i));
    end
    if b_correct==1
        binned_spks(bin_start_times==0)=0;
    end
    spk_density_fxn=conv(binned_spks,gaus_ker,'same');
    trl_fr(n,:)=spk_density_fxn./bin_size;
    
   if b_correct==0
       baseline(n)=mean(trl_fr(n,bin_start_times<-5));
       trl_fr(n,:)=trl_fr(n,:) - baseline(n);                 % possibly change from subtraction to division
   elseif b_correct==2
       baseline(n)=mean(trl_fr(n,bin_start_times<-5));
       trl_fr(n,:)=(trl_fr(n,:)-baseline(n))/baseline(n);
   else
       baseline(n)=NaN;
   end
end
mean_trl_fr=mean([trl_fr; trl_fr]);
std_fr=std(mean_trl_fr)./sqrt(size(trl_fr,1));

spk_d=plot(bin_start_times, mean_trl_fr, 'b', 'LineWidth',2.5); hold on
%plot(bin_start_times, 0, 'k--')
    
    