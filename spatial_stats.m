%
%
%        CALCULATES BACKGROUND STATE AND BACKGROUND COVARIANCE MATRIX
%
%
%
%
%==========================================================================
% DESCRIPTION
%==========================================================================
%
% The Background state and covariance matrices are calculated through
% multivariate statistics. To do so we suppose that decorrelated snapshots
% of mantle convection solution are stored in a folder. This script analyses
% these snapshots and calculates the mean state and covariance of the different
% components of the state (temperature field, velocity field...).
%
% The script is divided in 2 stages:
%     - Calculates the mean state, standard deviation and build the matrix
%     containing the successive states.
%     - calculates the correlations.
% The variants of the output are detailed in the section 'options'

clear all
close all

%==========================================================================
% OPTIONS
%==========================================================================

ti_start=1;
ti_end=404;

% On which fields does the correlations need to be performed?
stats_T=logical(1);
stats_Vs=logical(1);
stats_qs=logical(1);
stats_ages=logical(0);

% outputs marginal pdf?
marginal_pdf=true;

topTval=0;
symmetry=true;

geometry='Sphe';% 'YY'

%==========================================================================
% PATHS (GENERIC NAMES)
%==========================================================================
addpath('functions/')
% path and name for the long run:
iname='plates2';
ipath=strcat('/home/marie/Research/project2/initial_stats/plates2/dec_states/');

% Generic path of the results:
[opath]= check_or_create('/home/marie/Research/project2/initial_stats/plates2/','spatial2',false);


%==========================================================================
% CALCULATE MEAN AND STANDARD DEVIATIONS
%==========================================================================


[byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,iti_ad,erupta_total,botT_val,th_coord,ph_coord,r_coord]=catch_header(ipath,iname,'t',ti_start);
ntime=ti_end-ti_start+1;
nel=nthtot*nphtot*nrtot*nblocks;
nobs=nthtot*nphtot*nblocks;

if stats_T
    BT=struct;
    Tvec=zeros(ntime,nel);
end
if stats_Vs
    Vsvec=zeros(ntime,nobs);
    BVphs=struct;
end
if stats_qs
    Bq=struct;
    qsvec=zeros(ntime,nobs);
end
if stats_ages
    Bage=struct;
    agevec=zeros(ntime,nobs);
end


for k=1:ntime
    if stats_T
        [T]=read_field(ipath,iname,'t',k,geometry);
        Tvec(k,:)=T(:)';
        [BT]=add_to_bg(BT,T,k);
    end
    if stats_Vs
        [Vth,Vph,Vr,P,~]=read_vp(ipath,iname,'vp',k,geometry);
        Vphs=Vph(:,:,end);
        Vsvec(k,:)=Vphs(:)';
        [BVphs]=add_to_bg(BVphs,Vphs,k);
    end
    if stats_qs
        [q]=t_to_q(T,1,rg,nrtot,topTval,'top');
        qsvec(k,:)=q(:)';
        [Bq]=add_to_bg(Bq,q,k);
    end
    if stats_ages
        [age]=read_field(ipath,iname,'age',k,geometry);
        ages=age(:,:,end);
        agevec(k,:)=ages(:)';
        [Bage]=add_to_bg(Bage,q,k);
    end
    
end

if stats_T
    BT.sd=(Tvec-repmat(BT.m(:),[1,ntime])')'*((Tvec-repmat(BT.m(:),[1,ntime])')/(ntime-1)).^0.5;
end
if stats_Vs
    BVphs.sd=(Vsvec-repmat(BVphs.m(:),[1,ntime])')'*((Vsvec-repmat(BVphs.m(:),[1,ntime])')/(ntime-1)).^0.5;
end
if stats_qs
    Bq.sd=(qsvec-repmat(Bq.m(:),[1,ntime])')'*((qsvec-repmat(Bq.m(:),[1,ntime])')/(ntime-1)).^0.5;
end
if stats_ages
    Bage.sd=(agevec-repmat(Bage.m(:),[1,ntime])')'*((agevec-repmat(Bage.m(:),[1,ntime])')/(ntime-1)).^0.5;
end


Fields=struct;
[X,Y,Z,Ph,Th,R,RedFlag]=read_grid(iname,ipath,k,'t',geometry);
if stats_T
    M=[(1:ntime)'  BT.dm' BT.dsd'];
    header='frame number,average difference of mean  temperature,average difference of standard deviation of  temperature';
    dlmwrite(strcat(opath,'cvgce_mean_std_devT.txt'),header,'delimiter','');
    dlmwrite(strcat(opath,'cvgce_mean_std_devT.txt'),M,'delimiter',',','-append');
    if symmetry
        [BT]=sym_bg(BT,ntime*nphtot);
    end
    Fields.MeanT=BT.m;
    Fields.SigmaT=BT.sd;
    if strcmp(geometry,'Sphe')
        [~]=mat_to_pw(opath,iname,1,Fields,X,Y,Z);
    end
    M=[ 1-r_coord squeeze(BT.m(1,1,:)) squeeze(BT.sd(1,1,:)) squeeze(BT.msq(1,1,:)) squeeze(BT.rms(1,1,:))];
    header='radial coordinates, Temperature average, Temperature standard deviation, Temperature mean square, Temperature root mean square';
    dlmwrite(strcat(opath,'Temperature_stat.txt'),header,'delimiter','');
    dlmwrite(strcat(opath,'Temperature_stat.txt'),M,'delimiter',',','-append');
    save(strcat(opath,'BT'),'BT');
    
end
if stats_Vs
    M=[(1:ntime)'  BVphs.dm' BVphs.dsd'];
    header='frame number,average difference of mean  latitudinal velocity,average difference of standard deviation of latitudinal velocity';
    dlmwrite(strcat(opath,'cvgce_mean_std_devVth.txt'),header,'delimiter','');
    dlmwrite(strcat(opath,'cvgce_mean_std_devVth.txt'),M,'delimiter',',','-append');
    if symmetry
        [BVphs]=sym_bg(BVphs,ntime*nphtot);
    end
    M=[ BVphs.m(1,1)  BVphs.sd(1,1)  BVphs.msq(1,1)  BVphs.rms(1,1)];
    header='surface longitudinal velocity average, surface longitudinal velocity standard deviation, surface longitudinal velocity mean square, surface longitudinal velocity root mean square';
    dlmwrite(strcat(opath,'surface_longitudinal_velocity_stat.txt'),header,'delimiter','');
    dlmwrite(strcat(opath,'surface_longitudinal_velocity_stat.txt'),M,'delimiter',',','-append');
    save(strcat(opath,'BVphs'),'BVphs');
    
end
if stats_qs
    M=[(1:ntime)'  Bq.dm' Bq.dsd'];
    header='frame number,average difference of mean heat flux,average difference of standard deviation of heat flux';
    dlmwrite(strcat(opath,'cvgce_mean_std_devq.txt'),header,'delimiter','');
    dlmwrite(strcat(opath,'cvgce_mean_std_devq.txt'),M,'delimiter',',','-append');
    if symmetry
        [Bq]=sym_bg(Bq,ntime*nphtot);
    end
    M=[ Bq.m(1,1) Bq.sd(1,1) Bq.msq(1,1) Bq.rms(1,1)];
    header='surface heat flux average, surface heat flux standard deviation, surface heat flux mean square, surface heat flux root mean square';
    dlmwrite(strcat(opath,'heat_flux_stat.txt'),header,'delimiter','');
    dlmwrite(strcat(opath,'heat_flux_stat.txt'),M,'delimiter',',','-append');
    save(strcat(opath,'Bq'),'Bq');
    
end
if stats_ages
    M=[(1:ntime)'  Bage.dm' Bage.dsd'];
    header='frame number,average difference of mean age,average difference of standard deviation of age';
    dlmwrite(strcat(path,'cvgce_mean_std_devage.txt'),header,'delimiter','');
    dlmwrite(strcat(path,'cvgce_mean_std_devage.txt'),M,'delimiter',',','-append');
    if symmetry
        [Bage]=sym_bg(Bage,ntime*nphtot);
    end
    M=[ Bage.m(1,1) Bage.sd(1,1) Bage.msq(1,1) Bage.rms(1,1)];
    header='surface age average, surface age standard deviation, surface age mean square, surface age root mean square';
    dlmwrite(strcat(opath,'heat_flux_stat.txt'),header,'delimiter','');
    dlmwrite(strcat(opath,'heat_flux_stat.txt'),M,'delimiter',',','-append');
    
    save(strcat(opath,'Bage'),'Bage');
end

%==========================================================================
% CALCULATE STANDARD SCORE AND PROBABILITY DISTRIBUTIONS
%==========================================================================

if stats_T
    Tvec=(Tvec-repmat(BT.m(:)',[ntime,1]))./repmat(BT.sd(:)',[ntime,1]);
    if marginal_pdf
        % define min and max from state vector for temperature
        minrange=min(Tvec(:));
        maxrange=max(Tvec(:));
        % define amplitudes of intervals in which bins are counted
        deltarange=(maxrange-minrange)/1000;
        % define intervals in which bins are counted
        binrange=(minrange:deltarange:maxrange)';
        hist_tot=[];
        for ir=1:nrtot                                                     % loop over each radius (radial symmetry)
            if symmetry
                M=Tvec(:,(ir-1)*nphtot+1:ir*nphtot);                           % select all the values taken by temperature at a specific radius over time
                bincounts=histc(M(:),binrange);                                % counts number of temperature values for each interval
                bincounts=bincounts/deltarange/sum(bincounts);                 % density histogram, so that the integration of the surface below histogram is 1
                for i_range=1:length(binrange)
                    bincounts(i_range)=sum(bincounts(max(i_range-60,1):min(i_range+60,length(binrange))))/121;
                end
                gauss = 1/sqrt(2*pi)*exp(-1/2*(binrange).^2);
                header='Temperature value, empirical pdf,theoretical pdf';
                dlmwrite(strcat(opath,'Temperature_hist_',num2str(ir),'.txt'),header,'delimiter','');
                dlmwrite(strcat(opath,'Temperature_hist_',num2str(ir),'.txt'),[binrange bincounts gauss],'delimiter',',','-append');
                hist_tot=[hist_tot;(1-r_coord(ir))*ones(size(binrange)) binrange bincounts gauss]; % creates a table
            else
                for ith=1:nthtot
                    for iph=1:nphtot
                        M=Tvec(:,(ir-1)*nphtot*nthtot+(iph-1)*nthtot+ith);                           % select all the values taken by temperature at a specific radius over time
                        bincounts=histc(M(:),binrange);                                % counts number of temperature values for each interval
                        bincounts=bincounts/deltarange/sum(bincounts);                 % density histogram, so that the integration of the surface below histogram is 1
                        for i_range=1:length(binrange)
                            bincounts(i_range)=sum(bincounts(max(i_range-60,1):min(i_range+60,length(binrange))))/121;
                        end
                        gauss = 1/sqrt(2*pi)*exp(-1/2*(binrange).^2);
                        header='Temperature value, empirical pdf,theoretical pdf';
                        dlmwrite(strcat(opath,'Temperature_hist_',num2str(ith),num2str(iph),num2str(ir),'.txt'),header,'delimiter','');
                        dlmwrite(strcat(opath,'Temperature_hist_',num2str(ith),num2str(iph),num2str(ir),'.txt'),[binrange bincounts gauss],'delimiter',',','-append');
                    end
                end
            end
        end
        if symmetry
            dlmwrite(strcat(opath,'Temperature_hist_tot.txt'),hist_tot,'delimiter',',');
        end
    end
    
end
if stats_Vs
    Vsvec=(Vsvec-repmat(BVphs.m(:)',[ntime,1]))./repmat(BVphs.sd(:)',[ntime,1]);
    if marginal_pdf
        % define min and max from state vector for surface velocity
        minrange=min(Vsvec(:));
        maxrange=max(Vsvec(:));
        % define amplitudes of intervals in which bins are counted
        deltarange=(maxrange-minrange)/100;
        % define intervals in which bins are counted
        binrange=(minrange:deltarange:maxrange)';
        if symmetry
            bincounts=histc(Vsvec(:),binrange);
            bincounts=bincounts/deltarange/sum(bincounts);                 % density histogram, so that the integration of the surface below histogram is 1
            gauss = 1/sqrt(2*pi)*exp(-1/2*(binrange).^2);
            header='Velocity value, empirical pdf,theoretical pdf';
            dlmwrite(strcat(opath,'Velocity_hist.txt'),header,'delimiter','');
            dlmwrite(strcat(opath,'Velocity_hist.txt'),[binrange bincounts gauss],'delimiter',',','-append');
        else
            for ith=1:nthtot
                for iph=1:nphtot
                    bincounts=histc(Vsvec(:,(iph-1)*nthtot+ith),binrange);
                    bincounts=bincounts/deltarange/sum(bincounts);                 % density histogram, so that the integration of the surface below histogram is 1
                    gauss = 1/sqrt(2*pi)*exp(-1/2*(binrange).^2);
                    header='Velocity value, empirical pdf,theoretical pdf';
                    dlmwrite(strcat(opath,'Velocity_hist',num2str(ith),num2str(iph),'.txt'),header,'delimiter','');
                    dlmwrite(strcat(opath,'Velocity_hist',num2str(ith),num2str(iph),'.txt'),[binrange bincounts gauss],'delimiter',',','-append');
                end
            end
        end
    end
end
if stats_qs
    qsvec=(qsvec-repmat(Bq.m(:)',[ntime,1]))./repmat(Bq.sd(:)',[ntime,1]);
    if marginal_pdf
        minrange=min(qsvec(:));
        maxrange=max(qsvec(:));
        deltarange=(maxrange-minrange)/100;
        binrange=(minrange:deltarange:maxrange)';
        if symmetry
            bincounts=histc(qsvec(:),binrange);
            bincounts=bincounts/deltarange/sum(bincounts);                 % density histogram, so that the integration of the surface below histogram is 1
            gauss = 1/sqrt(2*pi)*exp(-1/2*(binrange).^2);
            header='Heat flux value, empirical pdf,theoretical pdf';
            dlmwrite(strcat(opath,'heat_flux_hist.txt'),header,'delimiter','');
            dlmwrite(strcat(opath,'heat_flux_hist.txt'),[binrange bincounts gauss],'delimiter',',','-append');
        else
            for ith=1:nthtot
                for iph=1:nphtot
                    bincounts=histc(qsvec(:, (iph-1)*nthtot+ith),binrange);
                    bincounts=bincounts/deltarange/sum(bincounts);                 % density histogram, so that the integration of the surface below histogram is 1
                    gauss = 1/sqrt(2*pi)*exp(-1/2*(binrange).^2);
                    header='Heat flux value, empirical pdf,theoretical pdf';
                    dlmwrite(strcat(opath,'heat_flux_hist',num2str(ith),num2str(iph),'.txt'),header,'delimiter','');
                    dlmwrite(strcat(opath,'heat_flux_hist',num2str(ith),num2str(iph),'.txt'),[binrange bincounts gauss],'delimiter',',','-append');
                    
                end
            end
        end
    end
end
if stats_ages
    agevec=(agevec-repmat(Bage.m(:)',[ntime,1]))./repmat(Bage.sd(:)',[ntime,1]);
    if marginal_pdf
        % Surface heat flux distribution
        dlmwrite(strcat(opath,'age_value.txt'),agevec(:),'delimiter',',');
        minrange=min(agevec(:));
        maxrange=max(agevec(:));
        deltarange=(maxrange-minrange)/100;
        binrange=(minrange:deltarange:maxrange)';
        if symmetry
            bincounts=histc(agevec(:),binrange);
            bincounts=bincounts/deltarange/sum(bincounts);                 % density histogram, so that the integration of the surface below histogram is 1
            gauss = 1/sqrt(2*pi)*exp(-1/2*(binrange).^2);
            header='age value, empirical pdf,theoretical pdf';
            dlmwrite(strcat(opath,'age_hist.txt'),header,'delimiter','');
            dlmwrite(strcat(opath,'age_hist.txt'),[binrange bincounts gauss],'delimiter',',','-append');
        else
            for ith=1:nthtot
                for iph=1:nphtot
                    bincounts=histc(agevec(:, (iph-1)*nthtot+ith),binrange);
                    bincounts=bincounts/deltarange/sum(bincounts);                 % density histogram, so that the integration of the surface below histogram is 1
                    gauss = 1/sqrt(2*pi)*exp(-1/2*(binrange).^2);
                    header='age value, empirical pdf,theoretical pdf';
                    dlmwrite(strcat(opath,'age_hist',num2str(ith),num2str(iph),'.txt'),header,'delimiter','');
                    dlmwrite(strcat(opath,'age_hist',num2str(ith),num2str(iph),'.txt'),[binrange bincounts gauss],'delimiter',',','-append');
                    
                end
            end
        end
    end
end


%==========================================================================
% CALCULATE CORRELATIONS
%==========================================================================


if stats_T
    if strcmp(geometry,'Sphe')
        
        if symmetry
            CorrTT=zeros(nthtot,nphtot,nrtot,nrtot);
            nsamp=zeros(nphtot,nrtot,nrtot);
            for robs=1:nrtot
                robs
                for i_dph=0:nphtot-1
                    for i_dr=0:nrtot-1
                        corr=0;
                        
                        for phobs=1:nphtot
                            if mod(phobs+i_dph,nphtot)==0
                                tphobs=nphtot;
                            else
                                tphobs=mod(phobs+i_dph,nphtot);
                            end
                            if mod(robs+i_dr,nrtot)==0
                                trobs=nrtot;
                            else
                                trobs=mod(robs+i_dr,nrtot);
                            end
                            nsamp(i_dph+1,robs,i_dr+1)=nsamp(i_dph+1,robs,i_dr+1)+ntime;
                            corr=corr+Tvec(:,(robs-1)*nphtot+phobs)'*Tvec(:,(trobs-1)*nphtot+tphobs);
                        end
                        CorrTT(1,i_dph+1,trobs,robs)=corr/(nsamp(i_dph+1,robs,i_dr+1)-1);
                    end
                end
            end
        else
            CorrTT=Tvec'*Tvec/(ntime-1);
            CorrTT=(CorrTT+CorrTT')/2;
        end
        save(strcat(opath,'CorrTT'),'CorrTT');
        FileNames= {};
        for ir=1:nrtot
            if symmetry
                CorrTTvis=CorrTT(1,:,:,ir);
                Fields=struct;
                Fields.Correlation=CorrTTvis;
                [fname_vtk]=mat_to_pw(opath,strcat(iname,'corrTT'),ir,Fields,X,Y,Z);
                FileNames{ir} = {r_coord(ir)+rcmb,1,fname_vtk};
            else
                for ith=1:nthtot
                    for iph=1:nphtot
                        CorrTTvis=reshape(CorrTT(:,(ir-1)*nphtot*nthtot+(iph-1)*nthtot+ith),[nthtot,nphtot,nrtot]);
                        Fields=struct;
                        Fields.Correlation=CorrTTvis;
                        [fname_vtk]=mat_to_pw(opath,strcat(iname,'corrTT'),(ir-1)*nphtot*nthtot+(iph-1)*nthtot+ith,Fields,X,Y,Z);
                        FileNames{(ir-1)*nphtot*nthtot+(iph-1)*nthtot+ith} = {r_coord(ir)+rcmb,1,fname_vtk};
                    end
                end
            end
        end
        
    end
    if stats_Vs
        if strcmp(geometry,'Sphe')
            if symmetry
                CorrTVphs=zeros(nthtot,nphtot,nrtot);
                for phobs=1:nphtot
                    for iph=1:phobs-1
                        for ir=1:nrtot
                            CorrTVphs(1,iph-phobs+nphtot+1,ir)=...
                                CorrTVphs(1,iph-phobs+nphtot+1,ir)+...
                                Vsvec(:,phobs)'*Tvec(:,(ir-1)*nphtot+iph);
                        end
                    end
                    for iph=phobs:nphtot
                        for ir=1:nrtot
                            CorrTVphs(1,iph-phobs+1,ir)=...
                                CorrTVphs(1,iph-phobs+1,ir)+...
                                Vsvec(:,phobs)'*Tvec(:,(ir-1)*nphtot+iph);
                        end
                    end
                end
                
                CorrTVphs=CorrTVphs/(ntime*nphtot-1);
            else
                CorrTVphs=Tvec'*Vsvec/(ntime-1);
            end
        end
        save(strcat(opath,'CorrTVphs'),'CorrTVphs');
        if symmetry
            Fields=struct;
            Fields.CorrelationTVphs=CorrTVphs;
            [~]=mat_to_pw(opath,strcat(iname,'corrTVphs'),1,Fields,X,Y,Z);
        else
            for ith=1:nthtot
                for iph=1:nphtot
                    Fields=struct;
                    Fields.CorrelationTVphs=reshape(CorrTVphs(:,(iph-1)*nthtot+ith),[nthtot,nphtot,nrtot]);
                    [~]=mat_to_pw(opath,strcat(iname,'corrTVphs'),(iph-1)*nthtot+ith,Fields,X,Y,Z);
                end
            end
        end
    end
    if stats_qs
        if strcmp(geometry,'Sphe')
            CorrTq=zeros(nthtot,nphtot,nrtot);
            for phobs=1:nphtot
                for iph=1:phobs-1
                    for ir=1:nrtot
                        CorrTq(1,iph-phobs+nphtot+1,ir)=...
                            CorrTq(1,iph-phobs+nphtot+1,ir)+...
                            qsvec(:,phobs)'*Tvec(:,(ir-1)*nphtot+iph);
                    end
                end
                for iph=phobs:nphtot
                    for ir=1:nrtot
                        CorrTq(1,iph-phobs+1,ir)=...
                            CorrTq(1,iph-phobs+1,ir)+...
                            qsvec(:,phobs)'*Tvec(:,(ir-1)*nphtot+iph);
                    end
                end
            end
            
            CorrTq=CorrTq/(ntime*nphtot-1);
            save(strcat(opath,'CorrTq'),'CorrTq');
            Fields=struct;
            Fields.CorrelationTq=CorrTq;
            [~]=mat_to_pw(opath,strcat(iname,'corrTq'),1,Fields,X,Y,Z);
            
        end
    end
    if stats_ages
        if strcmp(geometry,'Sphe')
            CorrTages=zeros(nthtot,nphtot,nrtot);
            for phobs=1:nphtot
                for iph=1:phobs-1
                    for ir=1:nrtot
                        CorrTages(1,iph-phobs+nphtot+1,ir)=...
                            CorrTages(1,iph-phobs+nphtot+1,ir)+...
                            agesvec(:,phobs)'*Tvec(:,(ir-1)*nphtot+iph);
                    end
                end
                for iph=phobs:nphtot
                    for ir=1:nrtot
                        CorrTages(1,iph-phobs+1,ir)=...
                            CorrTages(1,iph-phobs+1,ir)+...
                            agesvec(:,phobs)'*Tvec(:,(ir-1)*nphtot+iph);
                    end
                end
            end
            
            CorrTages=CorrTages/(ntime*nphtot-1);
            save(strcat(opath,'CorrTages'),'CorrTages');
            Fields=struct;
            Fields.CorrelationTages=CorrTages;
            [~]=mat_to_pw(opath,strcat(iname,'corrTages'),1,Fields,X,Y,Z);
        end
    end
end

if stats_Vs
    if strcmp(geometry,'Sphe')
        CorrVphsVphs=zeros(nthtot,nphtot);
        for phobs=1:nphtot
            for iph=1:phobs-1
                CorrVphsVphs(1,iph-phobs+nphtot+1)=...
                    CorrVphsVphs(1,iph-phobs+nphtot+1)+...
                    Vsvec(:,phobs)'*Vsvec(:,iph);
            end
            for iph=phobs:nphtot
                CorrVphsVphs(1,iph-phobs+1)=...
                    CorrVphsVphs(1,iph-phobs+1)+...
                    Vsvec(:,phobs)'*Vsvec(:,iph);
            end
        end
        
        CorrVphsVphs=CorrVphsVphs/(ntime*nphtot-1);
        save(strcat(opath,'CorrVphsVphs'),'CorrVphsVphs');
    end
    
    if stats_qs
        if strcmp(geometry,'Sphe')
            CorrqVphs=zeros(nthtot,nphtot);
            for phobs=1:nphtot
                for iph=1:phobs-1
                    CorrqVphs(1,iph-phobs+nphtot+1)=...
                        CorrqVphs(1,iph-phobs+nphtot+1)+...
                        qsvec(:,phobs)'*Vsvec(:,iph);
                end
                for iph=phobs:nphtot
                    CorrqVphs(1,iph-phobs+1)=...
                        CorrqVphs(1,iph-phobs+1)+...
                        qsvec(:,phobs)'*Vsvec(:,iph);
                end
            end
            
            CorrqVphs=CorrqVphs/(ntime*nphtot-1);
            save(strcat(opath,'CorrqVphs'),'CorrqVphs');
        end
    end
    if stats_ages
        if strcmp(geometry,'Sphe')
            CorragesVphs=zeros(nthtot,nphtot);
            for phobs=1:nphtot
                for iph=1:phobs-1
                    CorragesVphs(1,iph-phobs+nphtot+1)=...
                        CorragesVphs(1,iph-phobs+nphtot+1)+...
                        agesvec(:,phobs)'*Vsvec(:,iph);
                end
                for iph=phobs:nphtot
                    CorragesVphs(1,iph-phobs+1)=...
                        CorragesVphs(1,iph-phobs+1)+...
                        agesvec(:,phobs)'*Vsvec(:,iph);
                end
            end
            
            CorragesVphs=CorragesVphs/(ntime*nphtot-1);
            save(strcat(opath,'CorragesVphs'),'CorragesVphs');
        end
    end
end

if stats_qs
    if strcmp(geometry,'Sphe')
        Corrqq=zeros(nthtot,nphtot);
        for phobs=1:nphtot
            for iph=1:phobs-1
                Corrqq(1,iph-phobs+nphtot+1)=...
                    Corrqq(1,iph-phobs+nphtot+1)+...
                    qsvec(:,phobs)'*qsvec(:,iph);
            end
            for iph=phobs:nphtot
                Corrqq(1,iph-phobs+1)=...
                    Corrqq(1,iph-phobs+1)+...
                    qsvec(:,phobs)'*qsvec(:,iph);
            end
        end
        
        Corrqq=Corrqq/(ntime*nphtot-1);
        save(strcat(opath,'Corrqq'),'Corrqq');
    end
    if stats_ages
        if strcmp(geometry,'Sphe')
            Corrqages=zeros(nthtot,nphtot);
            for phobs=1:nphtot
                for iph=1:phobs-1
                    Corrqages(1,iph-phobs+nphtot+1)=...
                        Corrqages(1,iph-phobs+nphtot+1)+...
                        agesvec(:,phobs)'*qsvec(:,iph);
                end
                for iph=phobs:nphtot
                    Corrqages(1,iph-phobs+1)=...
                        Corrqages(1,iph-phobs+1)+...
                        agesvec(:,phobs)'*qsvec(:,iph);
                end
            end
            
            Corragesq=Corragesq/(ntime*nphtot-1);
            save(strcat(opath,'Corrqq'),'Corrqq');
        end
    end
end


