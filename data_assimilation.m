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

%==========================================================================
% PATHS (GENERIC NAMES)
%==========================================================================
addpath('functions/')
% path and name for the long run:
addpath('statistics/inst/')

iname='plates2';

% path of the observations
obspath='/home/marie/Research/project2/evolutions/plates2/evo2obs/qVobsnum3000dt10err0.01t1000/';
statpath='/home/marie/Research/project2/initial_stats/plates2/spatial/';
kpath=strcat('/home/marie/Research/project2/initial_stats/plates2/kalman',cnum(1/res),'/');
opath='/home/marie/Research/project1/data_assim/plates2/essai1/';
opathstag='/home/marie/Research/project1/data_assim/plates2/essai1/stagout/';
STAGex='../../../../../STAG/STAG_2013_02_08_marie/stag_dt/stagyympi';
np=8;
epath='/home/marie/Research/project2/initial_stats/plates2/dec_states/';
ename=iname;
enum=1;

%==========================================================================
% READ OBSERVATION DATABASE TIMEVEC
%==========================================================================

[~,timevec]=read_obs_header(obspath,'head')
ti_fn_vec=(timevec-1)/nwrite+1;

%==========================================================================
% ANALYSIS
%==========================================================================

load(strcat(statpath,'BT.mat'));
load(strcat(statpath,'Bq.mat'));
load(strcat(statpath,'BVphs.mat'));

nel=numel(BT.m);
nobs=numel(Bq.m);
[nthtot,nphtot,nrtot]=size(BT.m);
Bq.sdbis=Bq.sd;
BT.sdbis=BT.sd;
BVphs.sdbis=BVphs.sd;

load(strcat(kpath,'K.mat'));
[byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,iti_ad,erupta_total,botT_val,th_coord,ph_coord,r_coord]=catch_header(epath,ename,'t',enum);
[X,Y,Z,Ph,Th,R,RedFlag]=read_grid(ename,epath,enum,'t',geometry);
for iti_fn=1:length(ti_fn_vec)
    iti_fn
    [nty,corr,~]=read_obs_time(obspath,'meta',iti_fn);
    [val,err,nx,ny,nz,len,typ]=read_obs(obspath,'obs',iti_fn,2);
    [q]=val';
    [val,err,nx,ny,nz,len,typ]=read_obs(obspath,'obs',iti_fn,1);
    Vphs=val';
    
    % READ FORECAST STATE AND FORECAST OBSERVATIONS
    
    if iti_fn==1
        Tf=BT.m;
        Vphsf=BVphs.m;
        Vphf=zeros(size(BT.m));
        Vthf=zeros(size(BT.m));
        Vrf=zeros(size(BT.m));
        Pf=zeros(size(BT.m));
        qf=Bq.m;
        Vphsf=BVphs.m;
    else
        [byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,iti_ad,erupta_total,botTval,th_coord,ph_coord,r_coord]=catch_header(opathstag,iname,'vp',ti_fn_vec(iti_fn));
        [Tf]=read_field(opathstag,iname,'t',ti_fn_vec(iti_fn),geometry);
        [Vthf,Vphf,Vrf,Pf,~]=read_vp(opathstag,iname,'vp',ti_fn_vec(iti_fn),geometry);
        Vphsf=Vphf(:,:,end);
        [qf]=t_to_q(Tf,1,rg,nrtot,Ttop,'top');
    end
    
    % UPDATE OF COVARIANCE MATRIX Pa(n-1) => Pf(n)
    
    if iti_fn~=1
        meanqfmqo=sum(qf(:)-q(:))/nobs;
        Vqf=sum((qf(:)-q(:)-repmat(meanqfmqo,[nobs,1])).^2)/(nobs-1)+res*Bq.msq(1,1);
        Vnf=(rg(2*nrtot+1)-rg(2*nrtot))^2*Vqf;
        Vnm1f=BT.sdbis(1,1,nrtot)^2;
        fac=Vnf/Vnm1f
        BT.sdbis=sqrt(fac)*BT.sdbis;
        meanVphsfmVphso=sum(Vphsf(:)-Vphs(:))/nobs;
        VVf=sum((Vphsf(:)-Vphs(:)-repmat(meanVphsfmVphso,[nobs,1])).^2)/(nobs-1)+res*BVphs.msq(1,1);
        VVnm1f=BVphs.sdbis(1,1)^2;
        fac=VVf/VVnm1f
        BVphs.sdbis=sqrt(fac)*BVphs.sdbis;
    end
    Ta=Tf(:)+(K(1:nel,:)*((([q Vphs]-[qf Vphsf])./[Bq.sd BVphs.sd])')).*BT.sdbis(:);
    Ta=reshape(Ta,[nthtot,nphtot,nrtot]);
    Vphsa=Vphsf(:)+(K(nel+1:end,:)*((([q Vphs]-[qf Vphsf])./[Bq.sd BVphs.sd])')).*BVphs.sdbis(:);
    Vphsa=reshape(Vphsa,[nthtot,nphtot]);
    
    
    % cuts the anomalies
    Ta(find(Ta>Tmax))=Tmax;
    Ta(find(Ta<Tmin))=Tmin;
    Vpha=Vphf;
    Vpha(:,:,end)=Vphsa;
    
    % CREATES STAG INPUTS
    if iti_fn==1
        [byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,iti_ad,erupta_total,botT_val,th_coord,ph_coord,r_coord]=catch_header(epath,ename,'t',enum);
        iti_step=timevec(iti_fn);
        iti_ad=0;
    else
        [byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,iti_ad,erupta_total,botT_val,th_coord,ph_coord,r_coord]=catch_header(opathstag,iname,'t',ti_fn_vec(iti_fn));
    end
    nnth=1;
    nnph=1;
    nnr=1;
    write_scalar_file(opathstag,iname,'t',ti_fn_vec(iti_fn),Ta,file_format,nmagic,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,iti_ad,erupta_total,botT_val,th_coord,ph_coord,r_coord)
    write_scalar_file(opathstag,strcat(iname,'for'),'t',ti_fn_vec(iti_fn),Tf,file_format,nmagic,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,iti_ad,erupta_total,botT_val,th_coord,ph_coord,r_coord)
    if ti_fn_vec(iti_fn)==1
        [byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,iti_ad,erupta_total,botT_val,th_coord,ph_coord,r_coord]=catch_header(epath,ename,'vp',enum);
        iti_step=timevec(iti_fn);
        iti_ad=0;
    else
        [byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,iti_ad,erupta_total,botT_val,th_coord,ph_coord,r_coord]=catch_header(opathstag,iname,'vp',ti_fn_vec(iti_fn));
    end
    nnth=1;
    nnph=1;
    nnr=1;
    scalefac=1;
    write_vp_file(opathstag,iname,ti_fn_vec(iti_fn),Vthf,Vpha,Vrf,Pf,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,iti_ad,erupta_total,botT_val,th_coord,ph_coord,r_coord,scalefac)
    write_vp_file(opathstag,strcat(iname,'for'),ti_fn_vec(iti_fn),Vthf,Vpha,Vrf,Pf,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,iti_ad,erupta_total,botT_val,th_coord,ph_coord,r_coord,scalefac)
    
    
    
    % RUN STAG TILL THE NEXT ANALYSIS
    
    prog=pwd;
    cd(opathstag);
    % Read par file into cell A
    fid = fopen(strcat(opathstag,'par'),'r');
    i = 1;
    tline = fgetl(fid);
    A{i} = tline;
    while ischar(tline)
        i = i+1;
        tline = fgetl(fid);
        A{i} = tline;
    end
    fclose(fid);
    
    i=0;
    found=false;
    while ~found
        i=i+1;
        k = findstr(A{i},'nsteps');
        if ~isempty(k)
            found=true;
        end
    end
    if iti_fn<length(ti_fn_vec)
        nsteps=timevec(iti_fn+1);
        A{i} = ['        nsteps       =',num2str(nsteps)];
        
        % Write cell A into new par file
        fid = fopen(strcat(opathstag,'par'), 'w');
        for i = 1:numel(A)
            fprintf(fid,'%s\n',  A{i});
        end
        fclose(fid);
        
        disp('running STAG till the next observations')
        eval(strcat('! mpirun -np',32,num2str(np),' ./',STAGex,'>> a.txt'));
        disp('finished running STAG')
        cd(prog);
    end
end
