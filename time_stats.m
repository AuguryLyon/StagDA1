%
%
%        Computes Temporal evolution and corrélations
%
%
%
%
%==========================================================================
% DESCRIPTION
%==========================================================================
%
% From a path with binaries of:
%    - either a lyapunov experiment or...
%    - a single evolution,
% gives:
%    - the average global evolution of temperature, bottom heat flux, top heat
%    flux....
%    - the statistical equilibrium
%    - the lyapunov time
%    - time auto correlation.

clear all
close all

%==========================================================================
% OPTIONS
%==========================================================================

evo_type='simple'; % 'lyap'
num_evo=1; % number of evolutions computed
geometry='Sphe'; %'YY'

iname='plates2';
fname=iname;
% starting frame of the long run
ti_start=1;
ti_end=53546;
% On which fields does the correlations need to be performed?
stats_T=logical(1);
stats_Vs=logical(1);
stats_qs=logical(1);
stats_qb=logical(1);
stats_ages=logical(0);
stats_lneta=logical(1);

topTval=0;
botTval=0.9;
AUTOCORR=true;
samp_autoc=1000;
num_dat=10;

%==========================================================================
% PATHS (GENERIC NAMES)
%==========================================================================

addpath('functions/')
ipath='/home/marie/Research/project2/evolutions/plates2/evo2/';
opath='/home/marie/Research/project2/evolutions/plates2/evo2res/';
mkdir(opath)

%==========================================================================
% CALCULATES GLOBAL EVOLUTION OF PARAMETERS
%==========================================================================
evo=struct('ti_step',[],'ti_ad',[],'ti_fn',[]);
evo.iti_start=0;
% Vectors initialization
if stats_T
    evo.Tmin=[];
    evo.Tmean=[];
    evo.Tmax=[];
end
if stats_Vs
    evo.Vrmssurfmin=[];
    evo.Vrmssurfmean=[];
    evo.Vrmssurfmax=[];
end
if stats_lneta
    evo.lnetamin=[];
    evo.lnetamean=[];
    evo.lnetamax=[];
end
if stats_ages
    evo.agemin=[];
    evo.agemean=[];
    evo.agemax=[];
end
if stats_qs
    evo.qsmin=[];
    evo.qsmean=[];
    evo.qsmax=[];
end
if stats_qb
    evo.qbmin=[];
    evo.qbmean=[];
    evo.qbmax=[];
end
evo.time_evo=[];

for i_evo=1:num_evo
    for iti_fn=ti_start:ti_end
        ok=1;
        if stats_T
            if strcmp(evo_type,'lyap')
            fullnamet       = [ipath,fname,cnum(i_evo),'_t',cnum(iti_fn)]                    % name of the file to read
            else
            fullnamet       = [ipath,fname,'_t',cnum(iti_fn)]                    % name of the file to read
            end
            ok=ok*(exist(fullnamet,'file')==2);
        end
        if stats_Vs
            if strcmp(evo_type,'lyap')
            fullnamevp       = [ipath,fname,cnum(i_evo),'_vp',cnum(iti_fn)] ;                   % name of the file to read
            else
            fullnamevp       = [ipath,fname,'_vp',cnum(iti_fn)];                   % name of the file to read
            end
            ok=ok*(exist(fullnamevp,'file')==2);
        end
        if stats_lneta
            if strcmp(evo_type,'lyap')
            fullnameeta       = [ipath,fname,cnum(i_evo),'_eta',cnum(iti_fn)]   ;                 % name of the file to read
            else
            fullnameeta       = [ipath,fname,'_eta',cnum(iti_fn)];                    % name of the file to read
            end
            ok=ok*(exist(fullnameeta,'file')==2);
        end
        if  stats_ages
            if strcmp(evo_type,'lyap')
            fullnameage       = [ipath,fname,cnum(i_evo),'_age',cnum(iti_fn)]   ;                 % name of the file to read
            else
            fullnameage      = [ipath,fname,'_age',cnum(iti_fn)];                    % name of the file to read
            end
            ok=ok*(exist(fullnameage,'file')==2);
        end
        if ok
            evo(i_evo).ti_fn=[evo(i_evo).ti_fn ; iti_fn];
            [byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,iti_ad,erupta_total,botTval,th_coord,ph_coord,r_coord]=catch_header(ipath,fname,'vp',iti_fn);
            evo(i_evo).ti_ad=[evo(i_evo).ti_ad; iti_ad];
            evo(i_evo).ti_step=[evo(i_evo).ti_step; iti_step];
            [Vol]=calc_vol(th_coord,ph_coord,rcmb,rg,nthtot,nphtot,nrtot,nblocks,geometry,'');
            Vols=Vol(:,:,nrtot);
            Volb=Vol(:,:,1);
            if stats_T
                [T]=read_field(ipath,fname,'t',iti_fn,geometry);
                evo(i_evo).Tmean=[evo(i_evo).Tmean ; sum(T(:).*Vol(:))/sum(Vol(:))];
                evo(i_evo).Tmin= [evo(i_evo).Tmin ; min(T(:))];
                evo(i_evo).Tmax= [evo(i_evo).Tmax ; max(T(:))];
                
            end
            if stats_lneta
                [eta]=read_field(ipath,fname,'eta',iti_fn,geometry);
                eta=log(eta);
                evo(i_evo).lnetamean=[evo(i_evo).lnetamean ; sum(eta(:).*Vol(:))/sum(Vol(:))];
                evo(i_evo).lnetamin= [evo(i_evo).lnetamin ; min(eta(:))];
                evo(i_evo).lnetamax= [evo(i_evo).lnetamax ; max(eta(:))];
                
            end
            if stats_ages
                [age]=read_field(ipath,fname,'age',iti_fn,geometry);
                age=age(:,:,end);
                Volsage=Vols(age>0);
                age=age(age>0);
                
                evo(i_evo).agemean=[evo(i_evo).agemean ; sum(age(:).*Volsage(:))/sum(Volsage(:))];
                evo(i_evo).agemin= [evo(i_evo).agemin ; min(age(:))];
                evo(i_evo).agemax= [evo(i_evo).agemax ; max(age(:))];
            end
            if stats_Vs
                [Vth,Vph,~,~,scalefac]=read_vp(ipath,fname,'vp',iti_fn,geometry);
                norm=Vth(:,:,nrtot).^2+Vph(:,:,nrtot).^2;
                evo(i_evo).Vrmssurfmin=[ evo(i_evo).Vrmssurfmin ;...
                    min(sqrt(norm(:))) ];                            % root mean square velocity
                evo(i_evo).Vrmssurfmean=[ evo(i_evo).Vrmssurfmean ;...
                    sqrt(sum(norm(:).*Vols(:))/sum(Vols(:))) ];                            % root mean square velocity
                evo(i_evo).Vrmssurfmax=[ evo(i_evo).Vrmssurfmax ;...
                    max(sqrt(norm(:))) ];                            % root mean square velocity
            end
            if stats_qs
                if ~stats_T
                    [T]=read_field(ipath,fname,'t',iti_fn,geometry);
                end
                [qs]=t_to_q(T,1,rg,nrtot,topTval,'top');
                evo(i_evo).qsmean=[evo(i_evo).qsmean ; sum(qs(:).*Vols(:))/sum(Vols(:))];
                evo(i_evo).qsmin= [evo(i_evo).qsmin ; min(qs(:))];
                evo(i_evo).qsmax= [evo(i_evo).qsmax ; max(qs(:))];
            end
            if stats_qb
                if ~stats_T
                    [T]=read_field(ipath,fname,'t',iti_fn,geometry);
                end
                [qb]=t_to_q(T,1,rg,1,botTval,'bot');
                evo(i_evo).qbmean=[evo(i_evo).qbmean ; sum(qb(:).*Volb(:))/sum(Volb(:))];
                evo(i_evo).qbmin= [evo(i_evo).qbmin ; min(qb(:))];
                evo(i_evo).qbmax= [evo(i_evo).qbmax ; max(qb(:))];
            end
            
        end
    end
    
    header='frame number, adim time, timestep';
    evo(i_evo).time_evo=[evo(i_evo).ti_fn  evo(i_evo).ti_ad evo(i_evo).ti_step];
    if stats_T
        header=strcat(header,',min T, mean T, max T');
        evo(i_evo).time_evo=[evo(i_evo).time_evo  evo(i_evo).Tmin evo(i_evo).Tmean evo(i_evo).Tmax];
    end
    if stats_Vs
        header=strcat(header,',min Vs, mean Vs, max Vs');
        evo(i_evo).time_evo=[evo(i_evo).time_evo  evo(i_evo).Vrmssurfmin evo(i_evo).Vrmssurfmean evo(i_evo).Vrmssurfmax];
    end
    if stats_lneta
        header=strcat(header,',min lneta, mean lneta, max lneta');
        evo(i_evo).time_evo=[evo(i_evo).time_evo evo(i_evo).lnetamin evo(i_evo).lnetamean evo(i_evo).lnetamax ];
    end
    if stats_ages
        header=strcat(header,',min age, mean age, max age');
        evo(i_evo).time_evo=[evo(i_evo).time_evo evo(i_evo).agemin evo(i_evo).agemean evo(i_evo).agemax];
    end
    if stats_qs
        header=strcat(header,',min qs, mean qs, max qs');
        evo(i_evo).time_evo=[evo(i_evo).time_evo evo(i_evo).qsmin evo(i_evo).qsmean evo(i_evo).qsmax];
    end
    if stats_qb
        header=strcat(header,',min qb, mean qb, max qb');
        evo(i_evo).time_evo=[evo(i_evo).time_evo evo(i_evo).qbmin evo(i_evo).qbmean evo(i_evo).qbmax];
    end
    dlmwrite(strcat(opath,'time_evo',num2str(i_evo),'.csv'),header,'delimiter','');
    dlmwrite(strcat(opath,'time_evo',num2str(i_evo),'.csv'),evo(i_evo).time_evo,'delimiter',',','-append');
end

%==========================================================================
% DETERMINATION OF STATISTICAL EQUILIBRIUM
%==========================================================================
for i_evo=1:num_evo
    if stats_T
        fit_Tmin=zeros((length(evo(i_evo).time_evo(:,1))-num_dat),6);
        fit_Tmax=zeros((length(evo(i_evo).time_evo(:,1))-num_dat),6);
        fit_Tmean=zeros((length(evo(i_evo).time_evo(:,1))-num_dat),6);
    end
    if stats_Vs
        fit_Vrmssurfmin=zeros((length(evo(i_evo).time_evo(:,1))-num_dat),6);
        fit_Vrmssurfmean=zeros((length(evo(i_evo).time_evo(:,1))-num_dat),6);
        fit_Vrmssurfmax=zeros((length(evo(i_evo).time_evo(:,1))-num_dat),6);
    end
    if stats_lneta
        fit_lnetamin=zeros((length(evo(i_evo).time_evo(:,1))-num_dat),6);
        fit_lnetamax=zeros((length(evo(i_evo).time_evo(:,1))-num_dat),6);
        fit_lnetamean=zeros((length(evo(i_evo).time_evo(:,1))-num_dat),6);
    end
    if stats_ages
        fit_agemin=zeros((length(evo(i_evo).time_evo(:,1))-num_dat),6);
        fit_agemax=zeros((length(evo(i_evo).time_evo(:,1))-num_dat),6);
        fit_agemean=zeros((length(evo(i_evo).time_evo(:,1))-num_dat),6);
    end
    if stats_qs
        fit_qsmin=zeros((length(evo(i_evo).time_evo(:,1))-num_dat),6);
        fit_qsmax=zeros((length(evo(i_evo).time_evo(:,1))-num_dat),6);
        fit_qsmean=zeros((length(evo(i_evo).time_evo(:,1))-num_dat),6);
    end
    if stats_qb
        fit_qsmin=zeros((length(evo(i_evo).time_evo(:,1))-num_dat),6);
        fit_qsmax=zeros((length(evo(i_evo).time_evo(:,1))-num_dat),6);
        fit_qsmean=zeros((length(evo(i_evo).time_evo(:,1))-num_dat),6);
    end
    for iti_sta=1:length(evo(i_evo).time_evo(:,1))-num_dat
        ind=4;
        if stats_T
            [fit_Tmin(iti_sta,1:3)]=regression(evo(i_evo).time_evo(:,2),evo(i_evo).time_evo(:,ind),iti_sta,iti_sta+num_dat,1);
            ind=ind+1;
            [fit_Tmean(iti_sta,1:3)]=regression(evo(i_evo).time_evo(:,2),evo(i_evo).time_evo(:,ind),iti_sta,iti_sta+num_dat,1);
            ind=ind+1;
            [fit_Tmax(iti_sta,1:3)]=regression(evo(i_evo).time_evo(:,2),evo(i_evo).time_evo(:,ind),iti_sta,iti_sta+num_dat,1)	;
            ind=ind+1;
        end
        if stats_Vs
            [fit_Vrmssurfmin(iti_sta,1:3)]=regression(evo(i_evo).time_evo(:,2),evo(i_evo).time_evo(:,ind),iti_sta,iti_sta+num_dat,1);
            ind=ind+1;
            [fit_Vrmssurfmean(iti_sta,1:3)]=regression(evo(i_evo).time_evo(:,2),evo(i_evo).time_evo(:,ind),iti_sta,iti_sta+num_dat,1);
            ind=ind+1;
            [fit_Vrmssurfmax(iti_sta,1:3)]=regression(evo(i_evo).time_evo(:,2),evo(i_evo).time_evo(:,ind),iti_sta,iti_sta+num_dat,1);
            ind=ind+1;
        end
        if stats_lneta
            [fit_lnetamin(iti_sta,1:3)]=regression(evo(i_evo).time_evo(:,2),evo(i_evo).time_evo(:,ind),iti_sta,iti_sta+num_dat,1);
            ind=ind+1;
            [fit_lnetamean(iti_sta,1:3)]=regression(evo(i_evo).time_evo(:,2),evo(i_evo).time_evo(:,ind),iti_sta,iti_sta+num_dat,1);
            ind=ind+1;
            [fit_lnetamax(iti_sta,1:3)]=regression(evo(i_evo).time_evo(:,2),evo(i_evo).time_evo(:,ind),iti_sta,iti_sta+num_dat,1)	;
            ind=ind+1;
        end
        if stats_ages
            [fit_agemin(iti_sta,1:3)]=regression(evo(i_evo).time_evo(:,2),evo(i_evo).time_evo(:,ind),iti_sta,iti_sta+num_dat,1);
            ind=ind+1;
            [fit_agemean(iti_sta,1:3)]=regression(evo(i_evo).time_evo(:,2),evo(i_evo).time_evo(:,ind),iti_sta,iti_sta+num_dat,1);
            ind=ind+1;
            [fit_agemax(iti_sta,1:3)]=regression(evo(i_evo).time_evo(:,2),evo(i_evo).time_evo(:,ind),iti_sta,iti_sta+num_dat,1)	;
            ind=ind+1;
        end
        if stats_qs
            [fit_qsmin(iti_sta,1:3)]=regression(evo(i_evo).time_evo(:,2),evo(i_evo).time_evo(:,ind),iti_sta,iti_sta+num_dat,1);
            ind=ind+1;
            [fit_qsmean(iti_sta,1:3)]=regression(evo(i_evo).time_evo(:,2),evo(i_evo).time_evo(:,ind),iti_sta,iti_sta+num_dat,1);
            ind=ind+1;
            [fit_qsmax(iti_sta,1:3)]=regression(evo(i_evo).time_evo(:,2),evo(i_evo).time_evo(:,ind),iti_sta,iti_sta+num_dat,1)	;
            ind=ind+1;
        end
        if stats_qb
            [fit_qbmin(iti_sta,1:3)]=regression(evo(i_evo).time_evo(:,2),evo(i_evo).time_evo(:,ind),iti_sta,iti_sta+num_dat,1);
            ind=ind+1;
            [fit_qbmean(iti_sta,1:3)]=regression(evo(i_evo).time_evo(:,2),evo(i_evo).time_evo(:,ind),iti_sta,iti_sta+num_dat,1);
            ind=ind+1;
            [fit_qbmax(iti_sta,1:3)]=regression(evo(i_evo).time_evo(:,2),evo(i_evo).time_evo(:,ind),iti_sta,iti_sta+num_dat,1)	;
            ind=ind+1;
        end
    end
    
    evo(i_evo).iti_start =ti_start;
    if stats_T
        stat_eqTmin=find(fit_Tmin(:,1)<fit_Tmin(1,1)/200,1,'first');
        stat_eqTmean=find(fit_Tmean(:,1)<fit_Tmean(1,1)/200,1,'first');
        stat_eqTmax=find(fit_Tmax(:,1)<fit_Tmax(1,1)/200,1,'first');
        evo(i_evo).iti_start=max([evo(i_evo).iti_start,stat_eqTmin,stat_eqTmean,stat_eqTmax]);
    end
    if stats_Vs
        stat_eqVrmssurfmin=find(fit_Vrmssurfmin(:,1)<fit_Vrmssurfmin(1,1)/200,1,'first');
        stat_eqVrmssurfmean=find(fit_Vrmssurfmean(:,1)<fit_Vrmssurfmean(1,1)/200,1,'first');
        stat_eqVrmssurfmax=find(fit_Vrmssurfmax(:,1)<fit_Vrmssurfmax(1,1)/200,1,'first');
        evo(i_evo).iti_start=max([evo(i_evo).iti_start,stat_eqVrmssurfmin,stat_eqVrmssurfmean,stat_eqVrmssurfmax]);
    end
    if stats_lneta
        stat_eqlnetamin=find(fit_lnetamin(:,1)<fit_lnetamin(1,1)/200,1,'first');
        stat_eqlnetamean=find(fit_lnetamean(:,1)<fit_lnetamean(1,1)/200,1,'first');
        stat_eqlnetamax=find(fit_lnetamax(:,1)<fit_lnetamax(1,1)/200,1,'first');
        evo(i_evo).iti_start=max([evo(i_evo).iti_start,stat_eqlnetamin,stat_eqlnetamean,stat_eqlnetamax]);
    end
    if stats_ages
        stat_eqagemin=find(fit_agemin(:,1)<fit_agemin(1,1)/200,1,'first');
        stat_eqagemean=find(fit_agemean(:,1)<fit_agemean(1,1)/200,1,'first');
        stat_eqagemax=find(fit_agemax(:,1)<fit_agemax(1,1)/200,1,'first');
        evo(i_evo).iti_start=max([evo(i_evo).iti_start,stat_eqagemin,stat_eqagemean,stat_eqagemax]);
    end
    if stats_qs
        stat_eqqsmin=find(fit_qsmin(:,1)<fit_qsmin(1,1)/200,1,'first');
        stat_eqqsmean=find(fit_qsmean(:,1)<fit_qsmean(1,1)/200,1,'first');
        stat_eqqsmax=find(fit_qsmax(:,1)<fit_qsmax(1,1)/200,1,'first');
        evo(i_evo).iti_start=max([evo(i_evo).iti_start,stat_eqqsmin,stat_eqqsmean,stat_eqqsmax]);
    end
    if stats_qb
        stat_eqqbmin=find(fit_qbmin(:,1)<fit_qbmin(1,1)/200,1,'first');
        stat_eqqbmean=find(fit_qbmean(:,1)<fit_qbmean(1,1)/200,1,'first');
        stat_eqqbmax=find(fit_qbmax(:,1)<fit_qbmax(1,1)/200,1,'first');
        evo(i_evo).iti_start=max([evo(i_evo).iti_start,stat_eqqbmin,stat_eqqbmean,stat_eqqbmax]);
    end    
    dlmwrite(strcat(opath,'stat_eq',num2str(i_evo),'.csv'),evo(i_evo).iti_start,'delimiter','');
end

%==========================================================================
% LYAPUNOV TIME
%==========================================================================

if strcmp(evo_type,'lyap')
    for i_evo=1:num_evo
        if i_evo==1
        ti_fn=evo(i_evo).ti_fn;
        ti_ad=evo(i_evo).ti_ad;
        else
            a=min(length(ti_fn),length(evo(i_evo).ti_fn));
            ti_fn=ti_fn(1:a);
            ti_ad=ti_ad(1,a);
        end
    end
    err=zeros(1,length(ti_fn));
    for iti=1:length(ti_fn)
        [byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,iti_ad,erupta_total,botTval,th_coord,ph_coord,r_coord]=catch_header(ipath,strcat(fname,cnum(i_evo)),'t',iti);
        [Vol]=calc_vol(th_coord,ph_coord,rcmb,rg,nthtot,nphtot,nrtot,nblocks,geometry,'');
        if strcmp(geometry,'Sphe')
            [X,Y,Z]=coord_to_cartgrid(th_coord,ph_coord,r_coord,rcmb);
        end
        Tmean=zeros(nthtot,nphtot*nblocks,nrtot);
        Tmeansq=zeros(nthtot,nphtot*nblocks,nrtot);
        for i_evo=1:num_evo
            [T]=read_field(ipath,strcat(fname,cnum(i_evo)),'t',iti,geometry);
            Tmean=Tmean+T;
            Tmeansq=Tmeansq+T.^2;
        end
        Tmean=Tmean/num_evo;
        Tmeansq=Tmeansq/num_evo;
        Tsd=sqrt(Tmeansq-Tmean.^2);
        Fields=struct;
        Fields.Tmean=Tmean;
        Fields.Tmeansq=Tmeansq;
        Fields.Tsd=Tsd;
        if strcmp(geometry,'Sphe')
            [~]=mat_to_pw(opath,iname,1,Fields,X,Y,Z);
        end
        err(iti)=sum(Tsd(:).*Vol(:))/sum(Vol(:).*Tmean(:));
    end
    dlmwrite(strcat(opath,'lyap.csv'),[ti_fn(:),ti_ad(:),err(:)],'delimiter','');
end
%==========================================================================
% AUTOCORRELATION
%==========================================================================

if AUTOCORR
    [byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,iti_ad,erupta_total,botTval,th_coord,ph_coord,r_coord]=catch_header(ipath,iname,'t',ti_start);
    
    if stats_T
        Auto_corrT=zeros(samp_autoc,nrtot);
        meansqT=zeros(1,nrtot);
        meanT=zeros(1,nrtot);
    end
    if stats_Vs
        Auto_corrV=zeros(samp_autoc,1);
        meansqV=0;
        meanV=0;
    end
    if stats_lneta
        Auto_corrlneta=zeros(samp_autoc,nrtot);
        meansqlneta=0;
        meanlneta=0;
    end
    if stats_ages
        Auto_corrages=zeros(samp_autoc,1);
        meansqages=0;
        meanages=0;
    end
    if stats_qs
        Auto_corrqs=zeros(samp_autoc,1);
        meansqqs=0;
        meanqs=0;
    end
    if stats_qb
        Auto_corrqb=zeros(samp_autoc,1);
        meansqqb=0;
        meanqb=0;
    end
    ok=1;
    % starting after statistical equilibrium
    ti_start=evo(i_evo).iti_start;
    counter=0;
    while ok
        k=ti_start+samp_autoc*counter+1;
        if stats_T
            if strcmp(evo_type,'lyap')
            fullnamet       = [ipath,fname,cnum(i_evo),'_t',cnum(k)]                    % name of the file to read
            else
            fullnamet       = [ipath,fname,'_t',cnum(k)]                    % name of the file to read
            end
            ok=ok*(exist(fullnamet,'file')==2);
        end
        if stats_Vs
            if strcmp(evo_type,'lyap')
            fullnamevp       = [ipath,fname,cnum(i_evo),'_vp',cnum(k)] ;                   % name of the file to read
            else
            fullnamevp       = [ipath,fname,'_vp',cnum(k)];                   % name of the file to read
            end
            ok=ok*(exist(fullnamevp,'file')==2);
        end
        if stats_lneta
            if strcmp(evo_type,'lyap')
            fullnameeta       = [ipath,fname,cnum(i_evo),'_eta',cnum(k)]   ;                 % name of the file to read
            else
            fullnameeta       = [ipath,fname,'_eta',cnum(k)];                    % name of the file to read
            end
            ok=ok*(exist(fullnameeta,'file')==2);
        end
        if  stats_ages
            if strcmp(evo_type,'lyap')
            fullnameage       = [ipath,fname,cnum(i_evo),'_age',cnum(k)]   ;                 % name of the file to read
            else
            fullnameage      = [ipath,fname,'_age',cnum(k)];                    % name of the file to read
            end
            ok=ok*(exist(fullnameage,'file')==2);
        end
        ok
        ok=ok*(k<=ti_end)
        if ok
            if stats_T
                [T]=read_field(ipath,iname,'t',k,geometry);
                meansqT=meansqT+reshape(sum(sum(T.^2,1),2),[1,nrtot]);
                meanT=meanT+reshape(sum(sum(T,1),2),[1,nrtot]);
            end
            if stats_Vs
                [Vth,Vph,Vr,P,~]=read_vp(ipath,iname,'vp',k,geometry);
                Vphs=Vph(:,:,end);
                meansqV=meansqV+sum(sum(Vphs.^2,1),2);
                meanV=meanV+sum(sum(Vphs,1),2);
            end
            if stats_lneta
                [eta]=read_field(ipath,iname,'eta',k,geometry);
                lneta=log(eta);
                meansqlneta=meansqlneta+reshape(sum(sum(lneta.^2,1),2),[1,nrtot]);
                meanlneta=meanlneta+reshape(sum(sum(lneta,1),2),[1,nrtot]);
            end
            if stats_ages
                [age]=read_field(ipath,iname,'age',k,geometry);
                meansqages=meansqages+sum(sum(age.^2,1),2);
                meanages=meanages+sum(sum(age,1),2);
            end
            if stats_qs
                [qs]=t_to_q(T,1,rg,nrtot,topTval,'top');
                meansqqs=meansqqs+sum(sum(qs.^2,1),2);
                meanqs=meanqs+sum(sum(qs,1),2);
            end
            if stats_qb
                [qb]=t_to_q(T,1,rg,1,botTval,'bot');
                meansqqb=meansqqb+sum(sum(qb.^2,1),2);
                meanqb=meanqb+sum(sum(qb,1),2);
            end
            counter=counter+1;
        end
    end
    if stats_T
        meansqT=meansqT/(nphtot*nthtot*nblocks*counter);
        meanT=meanT/(nphtot*nthtot*nblocks*counter);
        VarT=meansqT-meanT.^2;
    end
    if stats_Vs
        meansqV=meansqV/(nphtot*nthtot*nblocks*counter);
        meanV=meanV/(nphtot*nthtot*nblocks*counter);
        VarV=meansqV-meanV.^2;
    end
    if stats_lneta
        meansqlneta=meansqlneta/(nphtot*nthtot*nblocks*counter);
        meanlneta=meanlneta/(nphtot*nthtot*nblocks*counter);
        Varlneta=meansqlneta-meanlneta.^2;
    end
    if stats_ages
        meansqages=meansqages/(nphtot*nthtot*nblocks*counter);
        meanages=meanages/(nphtot*nthtot*nblocks*counter);
        Varages=meansqages-meanages.^2;
    end
    if stats_qs
        meansqqs=meansqqs/(nphtot*nthtot*nblocks*counter);
        meanqs=meanqs/(nphtot*nthtot*nblocks*counter);
        Varqs=meansqqs-meanqs.^2;
    end
    if stats_qb
        meansqqb=meansqqb/(nphtot*nthtot*nblocks*counter);
        meanqb=meanqb/(nphtot*nthtot*nblocks*counter);
        Varqb=meansqqb-meanqb.^2;
    end
    counter=counter-1;
    ti_ad=zeros(samp_autoc,1);
    for ite=1:counter
        for a=1:samp_autoc
            k=a+(ite-1)*samp_autoc+ti_start;
            [~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,ad,~,~,~,~,~]=catch_header(ipath,iname,'t',k);
                if a==1
                    ad0=ad
                end
            ti_ad(a)=ti_ad(a)+ad-ad0;
            if stats_T
            [T]=read_field(ipath,iname,'t',k,geometry);
                if a==1
                    T0=T;
                    dT0=T0-repmat(reshape(meanT,[1,1,nrtot]),[nthtot,nphtot*nblocks]);
                end
                dT=T-repmat(reshape(meanT,[1,1,nrtot]),[nthtot,nphtot*nblocks]);
                Auto_corrT(a,:)=Auto_corrT(a,:)+reshape(sum(sum(dT0.*dT,1),2),[1,nrtot]);
            end
            if stats_Vs
                [Vth,Vph,Vr,P,~]=read_vp(ipath,iname,'vp',k,geometry);
                Vphs=Vph(:,:,end);
                if a==1
                    Vphs0=Vphs;
                    dVphs0=Vphs0-repmat(meanV,[nthtot,nphtot*nblocks]);
                end
                dVphs=Vphs-repmat(meanV,[nthtot,nphtot*nblocks]);
                Auto_corrV(a)=Auto_corrV(a)+sum(sum(dVphs0.*dVphs));
            end
            if stats_lneta
            [eta]=read_field(ipath,iname,'eta',k,geometry);
                lneta=log(eta);
                if a==1
                    lneta0=lneta;
                    dlneta0=lneta0-repmat(reshape(meanlneta,[1,1,nrtot]),[nthtot,nphtot*nblocks]);
                end
                dlneta=lneta-repmat(reshape(meanlneta,[1,1,nrtot]),[nthtot,nphtot*nblocks]);
                Auto_corrlneta(a,:)=Auto_corrlneta(a,:)+reshape(sum(sum(dlneta0.*dlneta,1),2),[1,nrtot]);
            end
            if stats_ages
            [age]=read_field(ipath,iname,'age',k,geometry);
                if a==1
                    age0=age;
                    dage0=age0-repmat(meanages,[nthtot,nphtot*nblocks]);
                end
                dage=age-repmat(meanages,[nthtot,nphtot*nblocks]);
                Auto_corrages(a)=Auto_corrages(a)+sum(sum(dage0.*dage));
            end
            if stats_qs
                [qs]=t_to_q(T,1,rg,nrtot,topTval,'top');
                if a==1
                    qs0=qs;
                    dqs0=qs0-repmat(meanqs,[nthtot,nphtot*nblocks]);
                end
                dqs=qs-repmat(meanqs,[nthtot,nphtot*nblocks]);
                Auto_corrqs(a)=Auto_corrqs(a)+sum(sum(dqs0.*dqs));
            end
            if stats_qb
                [qb]=t_to_q(T,1,rg,1,botTval,'bot');
                if a==1
                    qb0=qb;
                    dqb0=qb0-repmat(meanqb,[nthtot,nphtot*nblocks]);
                end
                dqb=qb-repmat(meanqs,[nthtot,nphtot*nblocks]);
                Auto_corrqb(a)=Auto_corrqb(a)+sum(sum(dqb0.*dqb));
            end
        end
    end
    
    ti_ad=ti_ad/counter;
    if stats_T
        Auto_corrT=Auto_corrT/(nphtot*nthtot*nblocks*counter-1)./repmat(VarT,size(Auto_corrT(:,1)));
        for k=1:nrtot
            M=[(1:samp_autoc)',ti_ad(:),Auto_corrT(:,k)];
            header='nframe,adim time, autocorr T';
            dlmwrite(strcat(opath,'autocorr_timeT',num2str(k),'.txt'),header,'delimiter','');
            dlmwrite(strcat(opath,'autocorr_timeT',num2str(k),'.txt'),M,'delimiter',',','-append');
        end
        a_coord=(1:samp_autoc)/samp_autoc;
        [A,B,C]  = meshgrid(r_coord,a_coord,th_coord);
        Fields=struct;
        Fields.autoT=Auto_corrT;
        [~]=mat_to_pwstraight(opath,'autocorr_timeT',1,Fields,A,B,C);
        [Vol]=calc_vol(th_coord,ph_coord,rcmb,rg,nthtot,nphtot,nrtot,nblocks,geometry,'');
        Auto_corr=sum(Auto_corrT.*repmat(squeeze(Vol(1,1,:))',size(Auto_corrT(:,1))),2)./sum(repmat(squeeze(Vol(1,1,:))',size(Auto_corrT(:,1))),2);
        M=[(1:samp_autoc)',ti_ad(:),Auto_corr(:)];
        header='frame number,adimensional time, autocorr';
        dlmwrite(strcat(opath,'autocorr_timeTglob.txt'),header,'delimiter','');
        dlmwrite(strcat(opath,'autocorr_timeTglob.txt'),M,'delimiter',',','-append');
    end
    if stats_Vs
        Auto_corrV=Auto_corrV/(nphtot*nthtot*nblocks*counter-1)./repmat(VarV,size(Auto_corrV(:,1)));
        M=[(1:samp_autoc)',ti_ad(:),Auto_corrV(:)];
        header='frame number, adim time, autocorr V';
        dlmwrite(strcat(opath,'autocorr_timeV.txt'),header,'delimiter','');
        dlmwrite(strcat(opath,'autocorr_timeV.txt'),M,'delimiter',',','-append');
    end
    if stats_lneta
        Auto_corrlneta=Auto_corrlneta/(nphtot*nthtot*nblocks*counter-1)./repmat(Varlneta,size(Auto_corrlneta(:,1)));
        for k=1:nrtot
            M=[(1:samp_autoc)',ti_ad(:),Auto_corrlneta(:,k)];
            header='nframe,adim time, autocorr eta';
            dlmwrite(strcat(opath,'autocorr_timelneta',num2str(k),'.txt'),header,'delimiter','');
            dlmwrite(strcat(opath,'autocorr_timelneta',num2str(k),'.txt'),M,'delimiter',',','-append');
        end
        a_coord=(1:samp_autoc)/samp_autoc;
        [A,B,C]  = meshgrid(r_coord,a_coord,th_coord);
        Fields=struct;
        Fields.autolneta=Auto_corrlneta;
        [~]=mat_to_pwstraight(opath,'autocorr_timelneta',1,Fields,A,B,C);
        [Vol]=calc_vol(th_coord,ph_coord,rcmb,rg,nthtot,nphtot,nrtot,nblocks,geometry,'');
        Auto_corr=sum(Auto_corrlneta.*repmat(squeeze(Vol(1,1,:))',size(Auto_corrlneta(:,1))),2)./sum(repmat(squeeze(Vol(1,1,:))',size(Auto_corrlneta(:,1))),2);
        M=[(1:samp_autoc)',ti_ad(:),Auto_corr(:)];
        header='frame num,adimensional time, autocorr';
        dlmwrite(strcat(opath,'autocorr_timelnetaglob.txt'),header,'delimiter','');
        dlmwrite(strcat(opath,'autocorr_timelnetaglob.txt'),M,'delimiter',',','-append');
    end
    if stats_ages
        Auto_corrages=Auto_corrages/(nphtot*nthtot*nblocks*counter-1)./repmat(Varages,size(Auto_corrages(:,1)));
        M=[(1:samp_autoc)',ti_ad(:),Auto_corrages(:)];
        header='frame number,adim time, autocorr ages';
        dlmwrite(strcat(opath,'autocorr_timeages.txt'),header,'delimiter','');
        dlmwrite(strcat(opath,'autocorr_timeages.txt'),M,'delimiter',',','-append');
    end
    if stats_qs
        Auto_corrqs=Auto_corrqs/(nphtot*nthtot*nblocks*counter-1)./repmat(Varqs,size(Auto_corrqs(:,1)));
        M=[(1:samp_autoc)',ti_ad(:),Auto_corrqs(:)];
        header='frame number,adim time, autocorr qs';
        dlmwrite(strcat(opath,'autocorr_timeqs.txt'),header,'delimiter','');
        dlmwrite(strcat(opath,'autocorr_timeqs.txt'),M,'delimiter',',','-append');
    end
    if stats_qb
        Auto_corrqb=Auto_corrqb/(nphtot*nthtot*nblocks*counter-1)./repmat(Varqb,size(Auto_corrqb(:,1)));
        M=[(1:samp_autoc)', ti_ad(:),Auto_corrqb(:)];
        header='frame number, adim time, autocorr qb';
        dlmwrite(strcat(opath,'autocorr_timeqb.txt'),header,'delimiter','');
        dlmwrite(strcat(opath,'autocorr_timeqb.txt'),M,'delimiter',',','-append');
    end
end

