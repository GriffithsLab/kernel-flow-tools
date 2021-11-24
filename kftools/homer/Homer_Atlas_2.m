% Pipeline Kernel Homer3 and AtlasViewer

% NOTE: Beforehand do setpaths in Homer3 folder and add Atlas-viewer to the path.

%Then go to folder with digpts.txt, probe.SD and your snirf moments file
function [dc,fwmodel] = Homer_Atlas(snirfFileName,dirname_atlas_viewer,currFolder)

[fwmodel, probe] = ForwardModel(dirname_atlas_viewer, currFolder);
snirf_data = MeanMoments(snirfFileName);
%load('C:\Users\zheng_wang\Downloads\Homer3-master\Homer3-master\Test_John_FT\homerOutput\atlasViewer.mat');
dc = HomerProcessing(snirf_data, probe);
SaveData(dc,fwmodel.Adot);
%SaveData(dc);
end


% 1) Get Mean of the moments

function res = MeanMoments(snirfFileName)

snirf = SnirfLoad(snirfFileName); %original file
snirf_cop = snirf;
timeSeries = snirf.data.dataTimeSeries;
measurementList = snirf.data.measurementList;

snirf_cop.data.dataTimeSeries = timeSeries(:,2:3:length(timeSeries(1,:)));
snirf_cop.data.measurementList = measurementList(2:3:length(measurementList(1,:)));
%SnirfSave('../Test_ft_Renee/test_renee_ft.snirf', snirf_cop);
res = snirf_cop;
end

% 2) Homer3 processing

function dc = HomerProcessing(res, probe)
%[pname, fname] = fileparts(filename);
%if isempty(pname)
%    pname = pwd;
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load acquisition file 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Acquired data (See SNIRF spec for list of acquired data parameters) 
%acquired = SnirfClass([pname, '/', fname, '.snirf']);

acquired = res;
% String together user functions to get from acquired data 
% to HRF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data_dnan = acquired;
%tmp = acquired.data.dataTimeSeries;
%tmp(isnan(tmp))= mean(mean();
%data_dnan.data.dataTimeSeries = tmp;
acquired.data.dataTimeSeries = hmrR_PreprocessIntensity_NAN(res.data);

dod   = hmrR_Intensity2OD(acquired.data);
dod_nan_rm = dod;
dod_nan_rm.dataTimeSeries = hmrR_PreprocessIntensity_NAN(dod);
dod   = hmrR_BandpassFilt(dod_nan_rm, 0.01, 0.50);

dod_nan_rm = dod;
dod_nan_rm.dataTimeSeries = hmrR_PreprocessIntensity_NAN(dod);
dc    = hmrR_OD2Conc_zw(dod, acquired.probe, probe, [1,1]);
end



function dc = hmrR_OD2Conc_zw( dod, probe, probe_true, ppf )

dc = DataClass().empty();

for ii=1:length(dod)
    dc(ii) = DataClass();
    
    Lambda = probe.GetWls();
    SrcPos_data = probe.GetSrcPos();
    SrcPos = probe_true.srcpos;
    DetPos_data = probe.GetDetPos();
    DetPos = probe_true.detpos;
    nWav   = length(Lambda);
    ml     = dod(ii).GetMeasList();
    y      = dod(ii).GetDataTimeSeries();
    
    if length(ppf) < nWav
        warning('Length of ppf does not match the number of wavelengths. Falling back to ppf=1 for all wavelengths.');
        ppf = ones(1, nWav);
    elseif length(ppf) > nWav
        d = length(ppf)-nWav;
        ppf(end-d+1:end) = [];
    end
    
    if ~isempty(find(ppf==1))
        ppf = ones(size(ppf));
    end
    
    nTpts = size(y,1);
    
    e = GetExtinctions(Lambda);
    e = e(:,1:2) / 10; % convert from /cm to /mm
    einv = inv( e'*e )*e';
    %size(einv)
    
    
    for idx=1:length(probe_true.ml)
        k = 3*(idx-1)+1;
        %idx1 = lst(idx);
        %idx2 = find( probe.ml(:,4)>1 & probe.ml(:,1)== probe.ml(idx1,1) & ml(:,2)==ml(idx1,2) );
        
        rho = norm(SrcPos(probe_true.ml(idx,1),:)-DetPos(probe_true.ml(idx,2),:));
        
        idx1 = find(ml(:,4) ==1 & ml(:,1) == probe_true.ml(idx,1) &  ml(:,2)== probe_true.ml(idx,2) );
        idx2 = find(ml(:,4) ==2 & ml(:,1) == probe_true.ml(idx,1) &  ml(:,2)== probe_true.ml(idx,2) );
        if isempty(idx1)
            y_mean = ones(nTpts,2);
            y_mean(:,1)=mean(y,2);
            y_mean(:,2)=mean(y,2);
          if ppf(1)~=1
              
            y2(:,k:k+1) = ( einv * (y_mean./(ones(nTpts,1)*rho*ppf))' )';
          else
            y2(:,k:k+1) = ( einv * (y_mean./(ones(nTpts,1)))' )';
          end   
        else
          if ppf(1)~=1
            y2(:,k:k+1) = ( einv * (y(:,[idx1(1) idx2(1)'])./(ones(nTpts,1)*rho*ppf))' )';
          else
            y2(:,k:k+1) = ( einv * (y(:,[idx1(1) idx2(1)'])./(ones(nTpts,1)))' )';
          end
        end
        y2(:,k+2) = y2(:,k) + y2(:,k+1);
        dc(ii).AddChannelHbO(probe_true.ml(idx,1), probe_true.ml(idx,2));
        dc(ii).AddChannelHbR(probe_true.ml(idx,1), probe_true.ml(idx,2));
        dc(ii).AddChannelHbT(probe_true.ml(idx,1), probe_true.ml(idx,2));
    end   
    dc(ii).SetDataTimeSeries(y2);
    dc(ii).SetTime(dod(ii).GetTime());
end
end



% 3) AtlasViewer
function [fwmodel, probe] = ForwardModel(dirname_atlas_viewer, currFolder)

% Registration from Digpoints (+register probe to surface)

dirnameAtlas = [dirname_atlas_viewer, '/Data/Colin/'];
searchPaths = dirnameAtlas;

headvol = initHeadvol();
headsurf = initHeadsurf();
pialsurf = initPialsurf();
refpts = initRefpts();
digpts = initDigpts();
probe = initProbe();
fwmodel = initFwmodel();
labelssurf = initLabelssurf();
imgrecon = initImgRecon();

headvol    = getHeadvol(headvol, searchPaths);
headsurf   = getHeadsurf(headsurf, searchPaths);
pialsurf   = getPialsurf(pialsurf, searchPaths);

headsurf.pathname = dirnameAtlas;
headvol.pathname = dirnameAtlas;

refpts     = getRefpts(refpts, headsurf.pathname);
labelssurf = getLabelssurf(labelssurf, headsurf.pathname);
fwmodel    = getFwmodel(fwmodel, currFolder, pialsurf, headsurf, headvol, probe);
imgrecon   = getImgRecon(imgrecon, currFolder, fwmodel, pialsurf, probe);
  


%pwd = dirnameSubj
digpts     = getDigpts(digpts, currFolder, refpts);

probe = getProbe(probe, currFolder, digpts, headsurf, refpts);


% First determine transformation to monte carlo space from volume 
% Generate transformation from head volume to digitized points space
[rp_atlas, rp_subj] = findCorrespondingRefpts(refpts, digpts);
headvol.T_2digpts = gen_xform_from_pts(rp_atlas, rp_subj);
headvol.imgOrig = headvol.img;

% Register headvol to digpts but first check fwmodel if it's volume 
% is already registered to digpts. if it is then set the headvol object 
% to the fwmodel's headvol and reuse it. 
if ~isregisteredFwmodel(fwmodel, headvol)
    
    [headvol.img, digpts.T_2mc] = ...
        xform_apply_vol_smooth(headvol.img, headvol.T_2digpts);
    
    headvol.T_2mc   = digpts.T_2mc * headvol.T_2digpts;
    headvol.center = xform_apply(headvol.center, headvol.T_2mc);
    headvol.orientation = digpts.orientation;
    
    % The MC space volume changed invalidating fwmodel meshes and 
    % vol to surface mesh. We need to recalculate all of this.
    
   % fwmodel  = resetFwmodel(fwmodel, headvol);
   % imgrecon = resetImgRecon(imgrecon);
    
else
    
    % Reusing MC space headvol from fwmodel. 
    headvol = fwmodel.headvol;
    
    % We know that headvol.T_2mc = digpts.T_2mc * headvol.T_2digpts.
    % Here we need to recover digpts.T_2mc. We can do this from 
    % headvol.T_2mc and headvol.T_2digpts with a little matrix algebra
    digpts.T_2mc = headvol.T_2mc / headvol.T_2digpts;
    
end

digpts.refpts.pos = xform_apply(digpts.refpts.pos, digpts.T_2mc);
digpts.pcpos      = xform_apply(digpts.pcpos, digpts.T_2mc);
digpts.srcpos     = xform_apply(digpts.srcpos, digpts.T_2mc);
digpts.detpos     = xform_apply(digpts.detpos, digpts.T_2mc);
digpts.optpos     = [digpts.srcpos; digpts.detpos];
digpts.center     = digpts.refpts.center;

% Copy digitized optodes to probe object
if ~isempty(probe.optpos_reg)
    probe.optpos_reg = xform_apply(probe.optpos_reg, digpts.T_2mc * probe.T_2digpts);
else
    probe.optpos_reg = xform_apply(probe.optpos, digpts.T_2mc * probe.T_2digpts);
end

% move head surface to monte carlo space 
headsurf.mesh.vertices   = xform_apply(headsurf.mesh.vertices, headvol.T_2mc);
headsurf.center          = xform_apply(headsurf.center, headvol.T_2mc);
headsurf.centerRotation  = xform_apply(headsurf.centerRotation, headvol.T_2mc);

% move pial surface to monte carlo space 
pialsurf.mesh.vertices   = xform_apply(pialsurf.mesh.vertices, headvol.T_2mc);
pialsurf.center          = xform_apply(pialsurf.center, headvol.T_2mc);

% move anatomical labels surface to monte carlo space 
labelssurf.mesh.vertices = xform_apply(labelssurf.mesh.vertices, headvol.T_2mc);
labelssurf.center        = xform_apply(labelssurf.center, headvol.T_2mc);

% move ref points to monte carlo space 
refpts = xform_apply_Refpts(refpts, headvol.T_2mc);

% The fwmodel meshes are inherited from pial and head surf at init time. Therefore they are 
% in original unregistered volume space, not MC space. Therefore we have to transform it to 
% MC space. 
fwmodel.mesh.vertices       = xform_apply(fwmodel.mesh.vertices, headvol.T_2mc);
fwmodel.mesh_scalp.vertices = xform_apply(fwmodel.mesh_scalp.vertices, headvol.T_2mc);
fwmodel.mesh_orig.vertices  = xform_apply(fwmodel.mesh_orig.vertices, headvol.T_2mc);
fwmodel.mesh_scalp_orig.vertices = xform_apply(fwmodel.mesh_scalp_orig.vertices, headvol.T_2mc);

% The imgrecon meshes are inherited from pial surf at init time. Therefore they are 
% in original unregistered volume space, not MC space. Therefore we have to transform it to 
% MC space. 
imgrecon.mesh.vertices       = xform_apply(imgrecon.mesh.vertices, headvol.T_2mc);
imgrecon.mesh_orig.vertices  = xform_apply(imgrecon.mesh_orig.vertices, headvol.T_2mc);

% Finish Registration
if ~headsurf.isempty(headsurf)
    headobj = headsurf;
else
    headobj = headvol;
end

refpts.eeg_system.selected = '10-5';
refpts = set_eeg_active_pts(refpts, [], false);
% Finish registration
if isempty(probe.al)
    
    % Register probe by simply pulling (or pushing) optodes toward surface
    % toward (or away from) center of head.
    method = 'digpts';
    probe = pullProbeToHeadsurf(probe, headobj);
    probe.hOptodesIdx = 1;
   
else
    
    % Register probe using springs based method
    if headvol.isempty(headvol)
        menu('Error registering probe using spring relaxation. Headvol object is empty','OK');
        return;
    end
    method = 'springs';
    probe = probeRegisterSpringsMethod(probe, headvol, refpts);
  
end


% View registered optodes on the head surface
%probe = viewProbe(probe, 'registered');

% Draw measurement list and save handle
%probe = findMeasMidPts(probe);

%fwmodel = updateGuiControls_AfterProbeRegistration(probe, fwmodel, imgrecon, labelssurf);

%probe.hOptodesIdx = 1; 
%probe = setProbeDisplay(probe, headsurf, method);
% Precalculate
T_vol2mc = headvol.T_2mc;
fwmodel = genSensitivityProfileFromFluenceProf(fwmodel, probe, T_vol2mc, currFolder);

end

function SaveData(dc, Adot)
%output = load('homerOutput/test_john.mat');
%load('fw/Adot.mat');
%load('homerOutput/atlasViewer.mat')
%data = load('trials_markers');
%times = output.output.dc.time;

dc_HbO=dc.dataTimeSeries(:,1:3:length(dc.dataTimeSeries));
dc_HbO(isnan(dc_HbO))=0;
save('timeseries_HbO.mat', 'dc_HbO')
dc_HbR=dc.dataTimeSeries(:,2:3:length(dc.dataTimeSeries));
dc_HbR(isnan(dc_HbR))=0;
dc_HbR = dc_HbR;
save('timeseries_HbR.mat', 'dc_HbR')

alpha = 1e-6;
tmp=size(Adot);
n_channels = tmp(1);
n_vert = tmp(2);
A= Adot(:,:,1);
tmp= sort(diag(A*A'));
loc_min = find(tmp>0.001);
lambda= tmp(loc_min(1));
tmp= A'/(A*A'+lambda*eye(n_channels));
save('transform_channel2vertice_HbO.mat', 'tmp')
A= Adot(:,:,2);
tmp= sort(eig(A*A'));%sort(diag(A*A'));
loc_min = find(tmp>0.001);
lambda= tmp(loc_min(1));
tmp= A'/(A*A'+lambda*eye(n_channels));
save('transform_channel2vertice_HbR.mat', 'tmp')

end
