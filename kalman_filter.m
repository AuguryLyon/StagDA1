%
%
%        CALCULATES KALMAN FILTER AND TESTS
%
%
%
%
%==========================================================================
% DESCRIPTION
%==========================================================================
%
%  The Kalman filter is described as K=(HP)'[H(HP)'+R]^-1
%  Where H is the observation matrix, P is the background covariance
%  matrix, and R is the observation covariance matrix. several variants of
%  calculation of K have been developped (see options)
% The script is divided in 3 stages:
%     - Load/build H,P,R
%     - Calculates the Kalman gain matrix with different methods
%     - test the first analysis on samples
% The variants of the output are detailed in the section 'options'

clear all
close all
rng('shuffle')
%==========================================================================
% OPTIONS
%==========================================================================

% define the state: Temperature field, Velocity field, viscosity field....

state = 'TVs'; % 'TVs'


load_H=false;
load_R=false;

res=0.01;  % relative error on observation (same for all the observations)
vp=100000;
geometry='Sphe';
method='';


%==========================================================================
% PATHS (GENERIC NAMES)
%==========================================================================
addpath('functions/')
addpath('statistics/inst/')
% name of the experiment
iname='plates2';

ipath='/home/marie/Research/project2/initial_stats/plates2/spatial/';
opath=strcat('/home/marie/Research/project2/initial_stats/plates2/kalman',cnum(1/res),'/');
mkdir(opath);
epath='/home/marie/Research/project2/initial_stats/plates2/dec_states/';
ename=iname;
enum=1;

%==========================================================================
% BUILD H AND R
%==========================================================================

load(strcat(ipath,'BT.mat'),'BT');
load(strcat(ipath,'Bq.mat'),'Bq');
load(strcat(ipath,'BVphs.mat'),'BVphs');
[byte_offset,file_format,nval,nmagic,xyp,nthtot,nphtot,nrtot,nblocks,aspect,nnth,nnph,nnr,nnb,rg,rcmb,iti_step,iti_ad,erupta_total,botT_val,th_coord,ph_coord,r_coord]=catch_header(epath,ename,'t',enum);
[X,Y,Z,Ph,Th,R,RedFlag]=read_grid(ename,epath,enum,'t',geometry);

nel=nthtot*nphtot*nblocks*nrtot;
nobs=nthtot*nphtot*nblocks;
%OBSERVATION MATRIX H

V_mean=BVphs.m;
V_sigma=BVphs.sd;
Vy_sq=BVphs.msq;

q_mean=Bq.m;
q_sigma=Bq.sd;
q_sq=Bq.msq;

diag1=ones(nobs,1).*q_sq'./q_sigma.^2';
diag2=ones(nobs,1).*Vy_sq'./V_sigma.^2';
diag3=[diag1;diag2];
H=spalloc(2*nobs,nel+nobs,2*nobs);

for iph=1:nphtot
    for ith=1:nthtot
        ith
        H((iph-1)*nthtot+ith,nthtot*nphtot*(nrtot-1)+(iph-1)*nthtot+ith)=1;
    end
end

for iph=1:nphtot
    for ith=1:nthtot
        ith
        H(nobs+(iph-1)*nthtot+ith,nel+(iph-1)*nthtot+ith)=1;
    end
end
save(strcat(opath,'H.mat'),'H');


nit=0;

%==========================================================================
% CALCULATES KALMAN
%==========================================================================

load(strcat(ipath,'Pdeclil',num2str(vp),'.mat'));
R=res*spdiags(diag3,0,nobs*2,nobs*2);
save(strcat(opath,'R'),'R');
invR=inv(R);
HRH=H'*invR*H;

vHRHv=V'*HRH*V;

invL=inv(L);
invA=inv(invL+vHRHv);
A1=V'*H'*invR;

invAA1=invA*A1;

K=V*invAA1;
clear HRH
clear vHRHv
clear invL
clear invA
clear A1
clear invAA1

save(strcat(opath,'K'),'K');

for phobs=1:nobs+1:2*nobs
    Fields=struct;
    
    Fields.KR=reshape(K(1:nel,phobs),[nthtot,nphtot,nrtot]);
    [~]=mat_to_pw(opath,'Kalman',phobs,Fields,X,Y,Z);
end

