%
%
%        Selects decorrelated states from an evolution
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
% Creates a directory with snapshots of decorrelated times.


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

% time of statistical equilibrium
ti_steq=99;
% end time of experiment
ti_end=29334;
% adimensional time between 2 decorrelated states
dti_ad=0.05;

% starting number for the decorrelated states collection
start_num=1;

% Which fields to copy to the new directory?
copy_T=true;
copy_Vs=true;
copy_ages=false;
copy_eta=true;


%==========================================================================
% PATHS (GENERIC NAMES)
%==========================================================================

addpath('functions/')
ipath='/home/marie/Research/project1/long_runs/plates2/stag_out/';
opath='/home/marie/Research/project2/initial_stats/plates2/dec_states/';
mkdir(opath)


%==========================================================================
% Copy files
%==========================================================================


it=start_num;
for i_evo=1:num_evo
    %    time_evo=dlmread(strcat(opath,'stat_eq',num2str(i_evo),'.csv'),',',0,0);
    for iti_fn=ti_steq:ti_end
        ok=1;
        if copy_T
            if strcmp(evo_type,'lyap')
                fullnamet       = [ipath,iname,cnum(i_evo),'_t',cnum(iti_fn)]                    % name of the file to read
            else
                fullnamet       = [ipath,iname,'_t',cnum(iti_fn)]                    % name of the file to read
            end
            ok=ok*exist(fullnamet,'file');
        end
        if copy_Vs
            if strcmp(evo_type,'lyap')
                fullnamevp       = [ipath,iname,cnum(i_evo),'_vp',cnum(iti_fn)] ;                   % name of the file to read
            else
                fullnamevp       = [ipath,iname,'_vp',cnum(iti_fn)];                   % name of the file to read
            end
            ok=ok*exist(fullnamevp,'file');
        end
        if copy_eta
            if strcmp(evo_type,'lyap')
                fullnameeta       = [ipath,iname,cnum(i_evo),'_eta',cnum(iti_fn)]   ;                 % name of the file to read
            else
                fullnameeta       = [ipath,iname,'_eta',cnum(iti_fn)];                    % name of the file to read
            end
            ok=ok*exist(fullnameeta,'file');
        end
        if  copy_ages
            if strcmp(evo_type,'lyap')
                fullnameage       = [ipath,iname,cnum(i_evo),'_age',cnum(iti_fn)]   ;                 % name of the file to read
            else
                fullnameage      = [ipath,iname,'_age',cnum(iti_fn)];                    % name of the file to read
            end
            ok=ok*exist(fullnameage,'file');
        end
        if ok
            [byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,ti_ad,erupta_total,botT_val,th_coord,ph_coord,r_coord]=catch_header(ipath,iname,'vp',iti_fn);
            if it==start_num
                ti_ad_old=ti_ad;
            end
            idelta_ti=ti_ad- ti_ad_old;
            if idelta_ti>dti_ad || it==start_num                           % the timelaps is long enough to have decorrelated states
                ti_ad_old=ti_ad;
                if copy_T
                    fullnamet2=[opath,fname,'_t',cnum(it)]
                    copyfile(fullnamet,fullnamet2)
                end
                if copy_Vs
                    fullnamevp2=[opath,fname,'_vp',cnum(it)];
                    copyfile(fullnamevp,fullnamevp2)
                end
                if copy_eta
                    fullnameeta2=[opath,fname,'_eta',cnum(it)];
                    copyfile(fullnameeta,fullnameeta2)
                end
                if copy_ages
                    fullnameage2=[opath,fname,'_age',cnum(it)];
                    copyfile(fullnameage,fullnameage2)
                end
                it=it+1;
            end
        end
    end
end