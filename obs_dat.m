%
%
%     CREATES an Observation database for an evolution
%
%
%
%

clear all
close all
rng('shuffle')

%==========================================================================
% OPTIONS
%==========================================================================
ti_start=3000;
% ending frame of the reference experiment
time=1000;

% dimensional time in Ma between 2 inversions
ti_delta= 10;

res=0.01;  % relative error on observation (same for all the observations)

topTval=0;
geometry='Sphe';

observe_Vs=true;
observe_qs=true;
observe_ages=false;

vrms=477.89;
qrms=17.747;

%==========================================================================
% PATHS (GENERIC NAMES)
%==========================================================================
addpath('functions/')
% path and name for the long run:
addpath('statistics/inst/')


refname='plates2';
iname='plates2';
% path of the reference stag outputs
refpath=strcat('/home/marie/Research/project2/evolutions/plates2/evo2/');
%path of the results
respath=strcat('/home/marie/Research/project2/evolutions/plates2/evo2obs/');

if observe_ages
    respath=strcat(respath,'age');
end
if observe_qs
    respath=strcat(respath,'q');
end
if observe_Vs
    respath=strcat(respath,'V');
end
respath=strcat(respath,'obsnum',num2str(ti_start),'dt',num2str(ti_delta),'err',num2str(res),'t',num2str(time),'/');
mkdir(respath)
resdetpath=strcat(respath,'fullevo/');
mkdir(resdetpath)

ad_to_di=70*vrms;

filen=[refpath,refname,'_t',cnum(ti_start)];

[byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,ti_step,ti_ad,erupta_total,botT_val,th_coord,ph_coord,r_coord]=catch_header(refpath,refname,'t',ti_start);

t_adV=ti_ad;
ti_di=ti_ad*ad_to_di;
t_diV=ti_di;
ti_fnV=0;
ti_stepV=ti_step;
iti=1;
ti_samp=iti-1;

while exist(filen,'file') && t_diV(iti)-t_diV(1)<time
    ti_fnV=[ti_fnV, iti-1];
    [byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,ti_step,ti_ad,erupta_total,botT_val,th_coord,ph_coord,r_coord]=catch_header(refpath,refname,'t',ti_fnV(iti)+ti_start-1);
    t_adV=[t_adV,ti_ad];
    ti_di=ti_ad*ad_to_di;
    t_diV=[t_diV,ti_di];
    ti_stepV=[ti_stepV,ti_step];
    if t_diV(iti)-t_diV(ti_samp(end)+1)>ti_delta
        ti_samp=[ti_samp, iti-1];
    end
    filen2=[resdetpath,iname,'_t',cnum(ti_fnV(iti))];
    filenvp2=[resdetpath,iname,'_vp',cnum(ti_fnV(iti))];
    filenvp=[refpath,refname,'_vp',cnum(iti+ti_start-1)];
    copyfile(filen,filen2);
    copyfile(filenvp,filenvp2);
    iti=iti+1;
    filen=[refpath,refname,'_t',cnum(iti+ti_start-1)]
end

ti_stepV=ti_stepV-ti_stepV(1)+1;
ntimes=length(ti_samp);

write_obs_header(respath,'head',nmagic,ntimes,ti_stepV(ti_samp+1),file_format)

nobs=nthtot*nphtot*nblocks;

R=eye(nobs*2);
diag1=ones(nobs,1)*vrms*sqrt(res)
diag2=ones(nobs,1)*qrms*sqrt(res)

for iti_fn=1:length(ti_samp)
    nty=2;
    R=1;
    a=nobs*2;
    write_obs_time(respath,'meta',nmagic,iti_fn,nty,R,a,file_format)
    
    % READ TRUE STATE AT THE TIME OF ANALYSIS
    [Vth,Vph,Vr,P,scalefac]=read_vp(refpath,refname,'vp',ti_samp(iti_fn)+ti_start,geometry);
    [T]=read_field(refpath,refname,'t',ti_samp(iti_fn)+ti_start,geometry);
    q=t_to_q(T,1,rg,nrtot,topTval,'top');
    qt=q;
    % true observations
    
    Vphs=Vph(:,:,end);
    Vphst=Vphs;
    
    % introduce noise in true observation
    dY=randn(nobs,1).*diag1;
    Vphs=Vphs+reshape(dY(1:nobs),[nthtot,nphtot*nblocks]);
    write_obs(respath,'obs',nmagic,iti_fn,1,squeeze(Vphs),diag1,ones(1,nphtot),1:nphtot,ones(1,nphtot)*nrtot,nobs,'VS',file_format)
    
    dY=randn(nobs,1).*diag2;
    q=q+reshape(dY(1:nobs),[nthtot,nphtot*nblocks]);
    write_obs(respath,'obs',nmagic,iti_fn,2,squeeze(q),diag2,ones(1,nphtot),1:nphtot,ones(1,nphtot)*nrtot,nobs,'HF',file_format)
    [byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,ti_step,ti_ad,erupta_total,botT_val,th_coord,ph_coord,r_coord]=catch_header(refpath,refname,'t',ti_samp(iti_fn)+ti_start)
    nnth=1;
    nnph=1;
    nnr=1;
    
    write_scalar_file(respath,strcat(iname,'_true'),'t',iti_fn,T,file_format,nmagic,nthtot,...
        nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,1,0,erupta_total,...
        botT_val,th_coord,ph_coord,r_coord);
    
    [byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,ti_step,ti_ad,erupta_total,botT_val,th_coord,ph_coord,r_coord]=catch_header(refpath,refname,'vp',ti_samp(iti_fn)+ti_start)
    nnth=1;
    nnph=1;
    nnr=1;
    write_vp_file(respath,strcat(iname,'_true'),iti_fn,Vth,Vph,Vr,P,file_format,nval,nmagic,xyp,nthtot,...
        nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,1,0,erupta_total,...
        botT_val,th_coord,ph_coord,r_coord,scalefac);
    
    dlmwrite(strcat(respath,'velocity_',num2str(iti_fn),'.txt'),Vphs(:),'delimiter',',');
    dlmwrite(strcat(respath,'velocitytrue_',num2str(iti_fn),'.txt'),Vphst(:),'delimiter',',');
    
    dlmwrite(strcat(respath,'heatflux_',num2str(iti_fn),'.txt'),q(:),'delimiter',',');
    dlmwrite(strcat(respath,'heatfluxtrue_',num2str(iti_fn),'.txt'),qt(:),'delimiter',',');
    if strcmp(geometry,'Sphe')
        [X,Y,Z,Ph,Th,R,RedFlag]=read_grid(refname,refpath,ti_samp(iti_fn)+ti_start,'t',geometry);
        
        Fields=struct;
        Fields.T=T;
        Fields.Vph=Vph;
        Fields.Vr=Vr;
        Fields.Vth=Vth;
        
        [fname_vtk]=mat_to_pw(respath,iname,iti_fn,Fields,X,Y,Z);
    elseif strcmp(geometry,'YY')
    end
end
