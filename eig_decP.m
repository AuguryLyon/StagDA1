%
%
%       Performs eigen decomposition of matrix P
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

% define the state: Temperature field, Velocity field, viscosity field....
state_T=true;
state_vp=true;
state_age=false;

% conditionning of the final matrix
vp=100000;

geometry='Sphe';
method='';
load_P=true;
Pmatname='P.mat';
eigdec=true;
Pdecname='Pdecomp.mat';
%==========================================================================
% PATHS (GENERIC NAMES)
%==========================================================================
addpath('functions/')
% name of the experiment
iname='plates2';
enum=1;
ipath=strcat('/home/marie/Research/project2/initial_stats/plates2/dec_states/');%/';
%iname='plates2';


statpath='/home/marie/Research/project2/initial_stats/plates2/spatial/';
opath=statpath;

%==========================================================================
% build P
%==========================================================================


[byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,iti_ad,erupta_total,botT_val,th_coord,ph_coord,r_coord]=catch_header(ipath,iname,'t',enum);
[Vol]=calc_vol(th_coord,ph_coord,rcmb,rg,nthtot,nphtot,nrtot,nblocks,geometry,method);
[X,Y,Z,Ph,Th,R,RedFlag]=read_grid(iname,ipath,enum,'t',geometry);

nel=0;

if state_T
    nelT=nthtot*nphtot*nrtot*nblocks;
    nel=nel+nelT;
end
if state_vp
    nelVphs=nthtot*nphtot*nblocks;
    nel=nel+nelVphs;
end
if state_age
    nelage=nthtot*nphtot*nblocks;
    nel=nel+nelage;
end

if load_P
    
if eigdec
    load(strcat(statpath,Pdecname));
    
else
    load(strcat(statpath,Pmatname));
    disp('starting eigen decomposition')
    [Vbig,Lbig]=eig(P);
    disp('sorting eigendecomposition')
    [~,Ac]=sort(diag(Lbig),1,'descend');
    Lbig=Lbig(Ac,Ac);
    Vbig=Vbig(:,Ac);
    disp('saving eigendecomposition')
    save(strcat(opath,'Pdecomp.mat'),'Vbig','-v7.3');
    save(strcat(opath,'Pdecomp.mat'),'Lbig','-append','-v7.3');
    clear P

end    
else
if eigdec
    load(strcat(opath,'Pdecomp.mat'));
    
else
    if state_T
        load(strcat(statpath,'CorrTT'),'CorrTT');
        CorrTTmat=zeros(nelT,nelT);
        
        for r1=1:nrtot
            r1
            for ph1=1:nphtot
                for ph2=1:ph1-1
                    for r2=1:nrtot
                        CorrTTmat((r1-1)*nphtot+ph1,(r2-1)*nphtot+ph2)=CorrTT(1,ph2-ph1+nphtot+1,r2,r1);
                    end
                end
                for ph2=ph1:nphtot
                    for r2=1:nrtot
                        CorrTTmat((r1-1)*nphtot+ph1,(r2-1)*nphtot+ph2)=CorrTT(1,ph2-ph1+1,r2,r1);
                    end
                end
            end
        end
        clear CorrTT
        CorrTT=(CorrTTmat+CorrTTmat')/2;
        clear CorrTTmat
        if state_vp
            load(strcat(statpath,'CorrTVphs'),'CorrTVphs');
            load(strcat(statpath,'CorrVphsVphs'),'CorrVphsVphs');
            CorrVphsVphsmat=zeros(nelVphs,nelVphs);
            CorrTVphsmat=zeros(nelT,nelVphs);
            
            for ph1=1:nphtot
                for ph2=1:ph1-1
                    CorrVphsVphsmat(ph1,ph2)=CorrVphsVphs(1,ph2-ph1+nphtot+1);
                end
                for ph2=ph1:nphtot
                    CorrVphsVphsmat(ph1,ph2)=CorrVphsVphs(1,ph2-ph1+1);
                end
            end
            
            for r1=1:nrtot
                r1
                for ph1=1:nphtot
                    for ph2=1:ph1-1
                        CorrTVphsmat((r1-1)*nphtot+ph2,ph1)=CorrTVphs(1,ph2-ph1+nphtot+1,r1);
                    end
                    for ph2=ph1:nphtot
                        CorrTVphsmat((r1-1)*nphtot+ph2,ph1)=CorrTVphs(1,ph2-ph1+1,r1);
                    end
                end
            end
            
            clear CorrVphsVphs
            clear CorrTVphs
            CorrVphsVphs=(CorrVphsVphsmat+CorrVphsVphsmat')/2;
            clear CorrVphsVphsmat
            CorrTVphs=CorrTVphsmat;
            clear CorrTVphsmat
        end
        
        if state_age
            load(strcat(statpath,'CorrTage'),'CorrTage');
            load(strcat(statpath,'Corrageage'),'Corrageage');
            Corrageagemat=zeros(nelage,nelage);
            CorrTagemat=zeros(nelT,nelage);
            
            for ph1=1:nphtot
                for ph2=1:ph1-1
                    Corrageagemat(ph1,ph2)=Corrageage(1,ph2-ph1+nphtot+1);
                end
                for ph2=ph1:nphtot
                    Corrageagemat(ph1,ph2)=Corrageage(1,ph2-ph1+1);
                end
            end
            
            for r1=1:nrtot
                r1
                for ph1=1:nphtot
                    for ph2=1:ph1-1
                        CorrTagemat((r1-1)*nphtot+ph2,ph1)=CorrTage(1,ph2-ph1+nphtot+1,r1);
                    end
                    for ph2=ph1:nphtot
                        CorrTagemat((r1-1)*nphtot+ph2,ph1)=CorrTage(1,ph2-ph1+1,r1);
                    end
                end
            end
            
            clear Corrageage
            clear CorrTage
            Corrageage=(Corrageagemat+Corrageagemat')/2;
            clear Corrageagemat
            CorrTVphs=CorrTagemat;
            clear CorrTVphsmat
        end
        
        P=zeros(nel,nel);
        P(1:nelT,1:nelT)=CorrTT;
        clear CorrTT
        
        if state_vp
            P(1:nelT,nelT+1:nelT+nelVphs)=CorrTVphs;
            P(nelT+1:nelT+nelVphs,1:nelT)=CorrTVphs';
            P(nelT+1:nelT+nelVphs,nelT+1:nelT+nelVphs)=CorrVphsVphs;
            clear CorrTVphs
            clear CorrVphsVphs
        end
        
    end
    
    disp('starting eigen decomposition')
    [Vbig,Lbig]=eig(P);
    disp('sorting eigendecomposition')
    [~,Ac]=sort(diag(Lbig),1,'descend');
    Lbig=Lbig(Ac,Ac);
    Vbig=Vbig(:,Ac);
    disp('saving eigendecomposition')
    save(strcat(opath,'Pdecomp.mat'),'Vbig','-v7.3');
    save(strcat(opath,'Pdecomp.mat'),'Lbig','-append','-v7.3');
    clear P
end
end
disp(' eigendecomposition finished')

[ind]=find(diag(Lbig)<(Lbig(1,1)/vp),1,'first');
V=zeros(nel,ind);
for it1=1:ind
    V(:,it1)=Vbig(:,it1);
    if state_T
        T=reshape(V(1:nelT,it1),[nthtot,nphtot*nblocks,nrtot]);
        Vph=zeros(nthtot,nphtot*nblocks,nrtot);
        Vph(:,:,end)=reshape(V(nelT+1:nelT+nelVphs,it1),[nthtot,nphtot*nblocks]);
        Fields=struct;
        Fields.T=T;
        Fields.Vph=Vph;
        oname=strcat(iname,'eigenvec');
        if strcmp(geometry,'YY')
            [fname_number,fname_vtk]=create_vtu(opath,oname,X(:,1:nphtot,:),Y(:,1:nphtot,:),Z(:,1:nphtot,:),X(:,nphtot+1:end,:),Y(:,nphtot+1:end,:),Z(:,nphtot+1:end,:),Fields,ti_fn);
        else
            [fname_vtk]=mat_to_pw(opath,oname,it1,Fields,X,Y,Z);
        end
    end
end
L=Lbig(1:ind,1:ind);

clear Vbig
clear Lbig
save(strcat(opath,'Pdeclil',num2str(vp),'.mat'),'V','-v7.3');
save(strcat(opath,'Pdeclil',num2str(vp),'.mat'),'L','-append','-v7.3');
Proot=V*(L.^0.5);
save(strcat(opath,'Prootlil',num2str(vp),'.mat'),'Proot','-v7.3');
