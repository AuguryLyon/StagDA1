%
%
%        PERFORMS ONE DATA ASSIMILATION FOR A REFERENCE EVOLUTION AND A SET
%        OF PARAMETERS
%
%
%
%
%==========================================================================
% DESCRIPTION
%==========================================================================
%
%  this code performs automatic sequential data assimilation
%
% - First analysis
% - Running the experiment till the next step when observations are
% available
% - perform data assimilation on this timestep
% - .....
% - till the final time has been reached
clear all
close all
rng('shuffle')
%==========================================================================
% OPTIONS
%==========================================================================

vp=100000;

Tmin=0;                                                   % minimum possible temperature for the experiment
Tmax=1;                                                   % maximum possible temperature for the experiment
Vrmssurfmin=0;                                            % minimum amplitude of surface velocity for the experiment
Vrmssurfmax=300000;                                            % maximum amplitude of surface velocity for the experiment

geometry='Sphe';
method='';
Ttop=0;
nwrite=20;
res=0.01;

vrms=477.89;
calc_errT=true;
calc_biasV=true;
calc_histq=true;

%==========================================================================
% PATHS (GENERIC NAMES)
%==========================================================================
addpath('functions/')
% path and name for the long run:
addpath('statistics/inst/')

truename='plates2';
estname='plates2';

% path of the observations
truepath='/home/marie/Research/project2/evolutions/plates2/evo2obs/qVobsnum3000dt10err0.01t1000/fullevo/';
obspath='/home/marie/Research/project2/evolutions/plates2/evo2obs/qVobsnum3000dt10err0.01t1000/';
estpath='/home/marie/Research/project1/data_assim/plates2/essai1/stagout/';

respath='/home/marie/Research/project1/data_assim/plates2/essai1/results/';
mkdir(respath);
%==========================================================================
% READ OBSERVATION DATABASE TIMEVEC
%==========================================================================

[~,timevec]=read_obs_header(obspath,'head')
ti_fn_vec=(timevec-1)/nwrite+1;
ad_to_di=70*vrms;
%==========================================================================
% TIMELOOP ON EVOLUTION
%==========================================================================
i_ana=1;
ntime=ti_fn_vec(end)-ti_fn_vec(1)+1;
nana=length(ti_fn_vec);
[byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,iti_ad,erupta_total,botTval,th_coord,ph_coord,r_coord]=catch_header(truepath,truename,'vp',1);

% Initialization of comparison true state/
errTa=zeros(nana,1);
err1DTa=zeros(nrtot,nana);

errTf_nana=zeros(nana,1);
err1DTf_nana=zeros(nrtot,nana);

errTf=zeros(ntime,1);
err1DTf=zeros(nrtot,ntime);

errqa=zeros(nana,1);
errqf=zeros(ntime,1);
errqf_nana=zeros(nana,1);

errVphsa=zeros(nana,1);
errVphsf_nana=zeros(nana,1);
errVphsf=zeros(ntime,1);

innovationmeanqf=zeros(1,nana);
innovationmeanVf=zeros(1,nana);

residualTfmean=zeros(nrtot,nana);
residualTfsigma=zeros(nrtot,nana);

residualqfmean=zeros(1,nana);
residualqfsigma=zeros(1,nana);

residualVphsfmean=zeros(1,nana);
residualVphsfsigma=zeros(1,nana);


residualTamean=zeros(nrtot,nana);
residualTasigma=zeros(nrtot,nana);

residualqamean=zeros(1,nana);
residualqasigma=zeros(1,nana);

residualVphsamean=zeros(1,nana);
residualVphsasigma=zeros(1,nana);


for iti_fn=ti_fn_vec(1):ti_fn_vec(end)
    
    
    [byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,iti_ad,erupta_total,botTval,th_coord,ph_coord,r_coord]=catch_header(truepath,truename,'vp',iti_fn);
    if iti_fn==ti_fn_vec(1)
        ti_ad0=iti_ad;
    end
    iti_di=ad_to_di*(iti_ad-ti_ad0);
    if iti_fn==ti_fn_vec(1)
        [Vol]=calc_vol(th_coord,ph_coord,rcmb,rg,nthtot,nphtot,nrtot,nblocks,geometry,'');
        Vols=Vol(:,:,end);
        [X,Y,Z,Ph,Th,R,RedFlag]=read_grid(truename,truepath,iti_fn,'t',geometry);
    end
    [Tt]=read_field(truepath,truename,'t',iti_fn,geometry);
    [Vtht,Vpht,Vrt,Pt,~]=read_vp(truepath,truename,'vp',iti_fn,geometry);
    Vphst=Vpht(:,:,end);
    [qt]=t_to_q(Tt,1,rg,nrtot,Ttop,'top');
    
    [byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,iti_ad,erupta_total,botTval,th_coord,ph_coord,r_coord]=catch_header(estpath,estname,'vp',iti_fn);
    [Tf]=read_field(estpath,estname,'t',iti_fn,geometry);
    [Vthf,Vphf,Vrf,Pf,~]=read_vp(estpath,estname,'vp',iti_fn,geometry);
    Vphsf=Vphf(:,:,end);
    [qf]=t_to_q(Tf,1,rg,nrtot,Ttop,'top');
    
    if iti_fn==ti_fn_vec(i_ana)
        Ta=Tf;
        Vtha=Vthf;
        Vpha=Vphf;
        Vra=Vrf;
        Pa=Pf;
        Vphsa=Vphsf;
        qa=qf;
        
        [byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,iti_ad,erupta_total,botTval,th_coord,ph_coord,r_coord]=catch_header(estpath,strcat(estname,'for'),'vp',iti_fn);
        [Tf]=read_field(estpath,strcat(estname,'for'),'t',iti_fn,geometry);
        [Vthf,Vphf,Vrf,Pf,~]=read_vp(estpath,strcat(estname,'for'),'vp',iti_fn,geometry);
        Vphsf=Vphf(:,:,end);
        [qf]=t_to_q(Tf,1,rg,nrtot,Ttop,'top');
        
        [nty,corr,~]=read_obs_time(obspath,'meta',i_ana);
        [val,err,nx,ny,nz,len,typ]=read_obs(obspath,'obs',i_ana,2);
        [qo]=val';
        [val,err,nx,ny,nz,len,typ]=read_obs(obspath,'obs',i_ana,1);
        Vphso=val';
        
        
        errTa(i_ana,1)=sqrt(sum((Ta(:)-Tt(:)).^2.*Vol(:))/sum(Tt(:).^2.*Vol(:)));
        for ir=1:nrtot
            err1DTa(ir,i_ana)=sqrt(sum(sum((Ta(:,:,ir)-Tt(:,:,ir)).^2.*Vol(:,:,ir),1),2)/sum(sum(Tt(:,:,ir).^2.*Vol(:,:,ir),1),2));
        end
        
        errTf_nana(i_ana,1)=sqrt(sum((Tf(:)-Tt(:)).^2.*Vol(:))/sum(Tt(:).^2.*Vol(:)));
        for ir=1:nrtot
            err1DTf_nana(ir,i_ana)=sqrt(sum(sum((Tf(:,:,ir)-Tt(:,:,ir)).^2.*Vol(:,:,ir),1),2)/sum(sum(Tt(:,:,ir).^2.*Vol(:,:,ir),1),2));
        end
        
        errqa(i_ana,1)=sqrt(sum((qa(:)-qt(:)).^2.*Vols(:))/sum(qt(:).^2.*Vols(:)));
        errqf_nana(i_ana,1)=sqrt(sum((qf(:)-qt(:)).^2.*Vols(:))/sum(qt(:).^2.*Vols(:)));
        
        errVphsa(i_ana,1)=sqrt(sum((Vphsa(:)-Vphst(:)).^2.*Vols(:))/sum(Vphst(:).^2.*Vols(:)));
        errVphsf_nana(i_ana,1)=sqrt(sum((Vphsf(:)-Vphst(:)).^2.*Vols(:))/sum(Vphst(:).^2.*Vols(:)));
        
        innovationq=qf-qo;
        innovationVphs=Vphsf-Vphso;
        innovationmeanqf(i_ana)=sum(innovationq(:).*Vols(:))/sum(Vols(:));
        innovationmeanVf(i_ana)=sum(innovationVphs(:).*Vols(:))/sum(Vols(:));
        
        for i=1:nrtot
            residualTfmean(ir,i_ana)=sum(sum(Tf(:,:,ir)-Tt(:,:,ir),1),2)/nthtot/nphtot/nblocks;
            residualTfsigma(ir,i_ana)=sqrt(sum(sum((Tf(:,:,ir)-Tt(:,:,ir)-repmat(residualTfmean(ir,i_ana),[nthtot,nphtot*nblocks])).^2,1),2)/(nthtot*nphtot*nblocks-1));
        end
        
        residualqfmean(1,i_ana)=sum(qf(:)-qt(:))/nthtot/nphtot/nblocks;
        residualqfsigma(1,i_ana)=sqrt(sum((qf(:)-qt(:)-repmat(residualqfmean(1,i_ana),size(qf(:)))).^2)/(nthtot*nphtot*nblocks-1));
        
        residualVphsfmean(1,i_ana)=sum(Vphsf(:)-Vphst(:))/nthtot/nphtot/nblocks;
        residualVphsfsigma(1,i_ana)=sqrt(sum((Vphsf(:)-Vphst(:)-repmat(residualVphsfmean(i_ana),size(Vphsf(:)))).^2)/(nthtot*nphtot*nblocks-1));
        
        for i=1:nrtot
            residualTamean(ir,i_ana)=sum(sum(Ta(:,:,ir)-Tt(:,:,ir),1),2)/nthtot/nphtot/nblocks;
            residualTasigma(ir,i_ana)=sqrt(sum(sum((Ta(:,:,ir)-Tt(:,:,ir)-repmat(residualTamean(ir,i_ana),[nthtot,nphtot*nblocks])).^2,1),2)/(nthtot*nphtot*nblocks-1));
        end
        
        residualqamean(1,i_ana)=sum(qa(:)-qt(:))/nthtot/nphtot/nblocks;
        residualqasigma(1,i_ana)=sqrt(sum((qa(:)-qt(:)-repmat(residualqamean(i_ana),size(qf(:)))).^2)/(nthtot*nphtot*nblocks-1));
        
        residualVphsamean(1,i_ana)=sum(Vphsa(:)-Vphst(:))/nthtot/nphtot/nblocks;
        residualVphsasigma(1,i_ana)=sqrt(sum((Vphsa(:)-Vphst(:)-repmat(residualVphsamean(i_ana),size(Vphsf(:)))).^2)/(nthtot*nphtot*nblocks-1));
        
        
        Fields=struct;
        Fields.Tt=Tt;
        Fields.Tf=Tf;
        Fields.Ta=Ta;

        Fields.Vpht=Vpht;
        Fields.Vrt=Vrt;
        Fields.Vtht=Vtht;

        Fields.Vphf=Vphf;
        Fields.Vrf=Vrf;
        Fields.Vthf=Vthf;
        
        Fields.Vpha=Vpha;
        Fields.Vra=Vra;
        Fields.Vtha=Vtha;

        [fname_vtk]=mat_to_pw(respath,strcat(truename,'ana'),i_ana,Fields,X,Y,Z);
        
        % print forecast error for analysis times
        if iti_fn==ti_fn_vec(1)
            header='ti_fn, ti_step, ti_ad, ti_di, errTf, errqf, errVf';
            dlmwrite(strcat(respath,'errTfana.csv'),header,'delimiter','')
        end
        M=[iti_fn iti_step iti_ad iti_di  errTf(iti_fn)  errqf(iti_fn)  errVphsf(iti_fn)];
        dlmwrite(strcat(respath,'errfana.csv'),M,'delimiter',',','-append')

        % print analyzed error for analysis times

        if iti_fn==ti_fn_vec(1)
            header='ti_fn, ti_step, ti_ad, ti_di, errTa, errqa, errVa';
            dlmwrite(strcat(respath,'erra.csv'),header,'delimiter','')
        end
        M=[iti_fn iti_step iti_ad iti_di  errTa(i_ana)  errqa(i_ana)  errVphsa(i_ana)];
        dlmwrite(strcat(respath,'erra.csv'),M,'delimiter',',','-append')
        
        % print innovation mean for analysis times

        if iti_fn==ti_fn_vec(1)
            header='ti_fn, ti_step, ti_ad, ti_di, innovationmeanqf, innovationmeanVf';
            dlmwrite(strcat(respath,'innovation.csv'),header,'delimiter','')
        end
        M=[iti_fn iti_step iti_ad iti_di  innovationmeanqf(i_ana)  innovationmeanVf(i_ana)];
        dlmwrite(strcat(respath,'innovation.csv'),M,'delimiter',',','-append')
        
        % print residual q mean for analysis times

        if iti_fn==ti_fn_vec(1)
            header='ti_fn, ti_step, ti_ad, ti_di, res mean qf, res sigma qf, resmean qa, res sigma qa';
            dlmwrite(strcat(respath,'residualq.csv'),header,'delimiter','')
        end
        M=[iti_fn iti_step iti_ad iti_di  residualqfmean(i_ana) residualqfsigma(i_ana) residualqamean(i_ana) residualqasigma(i_ana)];
        dlmwrite(strcat(respath,'residualq.csv'),M,'delimiter',',','-append')

        % print residual Vphs mean for analysis times

        if iti_fn==ti_fn_vec(1)
            header='ti_fn, ti_step, ti_ad, ti_di, res mean Vf, res sigma Vf, resmean Va, res sigma Va';
            dlmwrite(strcat(respath,'residualV.csv'),header,'delimiter','')
        end
        M=[iti_fn iti_step iti_ad iti_di  residualVphsfmean(i_ana) residualVphsfsigma(i_ana) residualVphsamean(i_ana) residualVphsasigma(i_ana)];
        dlmwrite(strcat(respath,'residualV.csv'),M,'delimiter',',','-append')

        % print errors 1D profiles
        header='r_coord, err1DTf, err1DTa';
        dlmwrite(strcat(respath,'errT1Dana',cnum(i_ana),'.csv'),header,'delimiter','')
        M=[r_coord  err1DTf_nana(:,i_ana) err1DTa(:,i_ana)];
        dlmwrite(strcat(respath,'errT1Dana',cnum(i_ana),'.csv'),M,'delimiter',',','-append')
   
        % print residual 1D profiles
        header='r_coord, res mean Tf, res sigma Tf, res mean Ta, res sigma Ta';
        dlmwrite(strcat(respath,'resT1Dana',cnum(i_ana),'.csv'),header,'delimiter','')
        M=[r_coord  residualTfmean(:,i_ana) residualTfsigma(:,i_ana) residualTamean(:,i_ana) residualTasigma(:,i_ana)];
        dlmwrite(strcat(respath,'resT1Dana',cnum(i_ana),'.csv'),M,'delimiter',',','-append')

        i_ana=i_ana+1
    end
    
    errTf(iti_fn,1)=sqrt(sum((Tf(:)-Tt(:)).^2.*Vol(:))/sum(Tt(:).^2.*Vol(:)));
    for ir=1:nrtot
        err1DTf(ir,iti_fn)=sqrt(sum(sum((Tf(:,:,ir)-Tt(:,:,ir)).^2.*Vol(:,:,ir),1),2)/sum(sum(Tt(:,:,ir).^2.*Vol(:,:,ir),1),2));
    end
    errqf(iti_fn,1)=sqrt(sum((qf(:)-qt(:)).^2.*Vols(:))/sum(qt(:).^2.*Vols(:)));
    errVphsf(iti_fn,1)=sqrt(sum((Vphsf(:)-Vphst(:)).^2.*Vols(:))/sum(Vphst(:).^2.*Vols(:)));
    
    Fields=struct;
    Fields.Tt=Tt;
    Fields.Tf=Tf;
    
    Fields.Vpht=Vpht;
    Fields.Vrt=Vrt;
    Fields.Vtht=Vtht;
    
    Fields.Vphf=Vphf;
    Fields.Vrf=Vrf;
    Fields.Vthf=Vthf;
    
    [fname_vtk]=mat_to_pw(respath,strcat(truename,'for'),iti_fn,Fields,X,Y,Z);

    % print forecast error for each time frame
    
    if iti_fn==ti_fn_vec(1)
        header='ti_fn, ti_step, ti_ad, ti_di, errTf, errqf, errVf';
        dlmwrite(strcat(respath,'errTf.csv'),header,'delimiter','')
    end
    M=[iti_fn iti_step iti_ad iti_di  errTf(iti_fn)  errqf(iti_fn)  errVphsf(iti_fn)];
    dlmwrite(strcat(respath,'errTf.csv'),M,'delimiter',',','-append')
    
    header='r_coord, err1DTf';
    dlmwrite(strcat(respath,'errTf1D',cnum(iti_fn),'.csv'),header,'delimiter','')
    M=[r_coord  err1DTf(:,iti_fn)];
    dlmwrite(strcat(respath,'errTf1D',cnum(iti_fn),'.csv'),M,'delimiter',',','-append')
    
end