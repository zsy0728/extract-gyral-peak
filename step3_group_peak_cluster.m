%step3 groupwise cluster
%% 5level peaks clusters
clear
clc
warning off

datadir = '/media/songyao/songyao/result/HCP_gyral_peak_MNINonLinear/groupwise_peaks/';
Surf = vtkSurfRead([datadir,'groupwise_peaks_count_smooth.vtk']);
face=Surf.Face;
vtx=Surf.Vtx;

label = Surf.Pdata{1,1}.val;
lvl_bg_thr = -209;
merge_thr = 7;
lvl_fg_thr = -45;
level=0-label;
disp('Initializing');

neighbor=cell(size(vtx,2),1);
for i=1:size(face,2)
    neighbor{face(1,i)}=[neighbor{face(1,i)};face(2:3,i)];
    neighbor{face(2,i)}=[neighbor{face(2,i)};face([1 3],i)];
    neighbor{face(3,i)}=[neighbor{face(3,i)};face(1:2,i)];
end
for i=1:length(neighbor)
    neighbor{i}=unique(neighbor{i});
end
map=coConnectedComponent(neighbor, level<=lvl_bg_thr);
o_surf.Vtx = vtx;
o_surf.Face = face-1;
o_surf.Pdata{1}.val = map;
o_surf.Pdata{1}.name = 'component';

prevThr=lvl_bg_thr;
for thr=lvl_bg_thr:1:lvl_fg_thr
    region=level>prevThr & level<=thr;
    %iteratively expand
    seed=find(map>0);
    conflict=cell(max(map(:)));
    while ~isempty(seed)
        mapclaim=[];
        for vid=seed'
            tmp=find(region(neighbor{vid}));
            if ~isempty(tmp)
                mapclaim=[mapclaim;[neighbor{vid}(tmp),ones(length(tmp),1)*map(vid)]];
            end
        end
        nextseed=[];
        while ~isempty(mapclaim)
            vid=mapclaim(1,1);
            region(vid)=0;
            rid=unique(mapclaim(mapclaim(:,1)==vid,2));
            if length(rid)<=1 %none confliction
                map(vid)=rid;
                nextseed=[nextseed; vid];
            else %confliction
                %chose the most similar one
                %can just use the first one to speed up
                mindis=1e16;
                minrid=1;
                for tmprid=rid'
                    tmpdis=(level(vid)-mean(level(map==tmprid)));
                    if tmpdis<mindis
                        mindis=tmpdis;
                        minrid=tmprid;
                    end
                end
                map(vid)=minrid;
            end
            mapclaim(mapclaim(:,1)==vid,:)=[];
        end
        seed=nextseed;
    end
    
    %check possible merge
    if thr==1
    end
    map=coComponentMerge(neighbor, map, level, merge_thr);
    %cnnected components for the rest regions
    if sum(region)>0
        remain=coConnectedComponent(neighbor, region);
        map(region)=remain(region)+max(map);
    end
    border=coComponentBorder(neighbor, map);
    prevThr=thr;
end

%% chongpai map bianhao
table = tabulate(map);
row = find(table(:,2)<=20);
for m=1:length(row)
    bianhao = table(row(m),1);
    map(map==bianhao) = 0;
end
tablenew = tabulate(map);
for n=1:length(tablenew(:,1))
    map(map==tablenew(n,1)) = n-1;
end
disp(['area : ',num2str(sum(map~=0))])
disp(['cluster number : ',num2str(max(map))])
o_surf.Vtx = vtx;
o_surf.Face = face-1;
o_surf.Pdata{1,1}.val = map;
o_surf.Pdata{1,1}.name = 'map';
vtkSurfWrite([datadir,'groupwise_peaks_cluster_(',num2str(lvl_fg_thr),')-(',num2str(merge_thr),')-(',num2str(lvl_bg_thr),').vtk'],o_surf);

% inflate surface
infsurf = vtkSurfRead('/media/songyao/songyao/data/HCP_s900_data/display_surface/100206.inflated_MSMAll.32k_fs_LR.vtk');
infsurf.Face = face-1;
infsurf.Pdata=[];
infsurf.Pdata{1,1}.val = map;
infsurf.Pdata{1,1}.name = 'map';
vtkSurfWrite([datadir,'groupwise_peaks_cluster_inf.vtk'],infsurf);

%% expand clusters
clc
clear
warning off

surfdir = '/media/songyao/songyao/result/HCP_gyral_peak_MNINonLinear/groupwise_peaks/';
input_surface_fname = [surfdir,'groupwise_peaks_cluster_(-45)-(7)-(-209).vtk'];
output_surface_fname = [surfdir,'groupwise_peaks_cluster_expand.vtk'];

label_smooth_flag = 1;
Surf = vtkSurfRead(input_surface_fname);
surf_c = ReadSurf_2('/media/songyao/songyao/data/HCP_s900_data/T1_surface/100206.white_MSMAll.32k_fs_LR_c.vtk.ply.vtk',{},1);
new_label = [];
for vtxID = 1:size(Surf.Vtx,2)
    [Neighbor_ID,Ring_ID] = Search_Neighbor_ID(surf_c,vtxID-1,label_smooth_flag,Inf);
    Neighbor_ID = Neighbor_ID + 1;
    tmp_label = Surf.Pdata{1,1}.val(Neighbor_ID);
    tmp_unique_label=unique(tmp_label);
    if tmp_unique_label(1) == 0
        tmp_unique_label(1) = [];
    end
    if ~isempty(tmp_unique_label)
        tmp_max_label = [];
        tmp_max_num = 0;
        for j = 1:size(tmp_unique_label,1)
            if size(find(tmp_label == tmp_unique_label(j)),1) > tmp_max_num
                tmp_max_num = size(find(tmp_label == tmp_unique_label(j)),1);
                tmp_max_label = tmp_unique_label(j);
            end
        end
        new_label(vtxID) = tmp_max_label;
    else
        new_label(vtxID) = 0;
    end
    
end
Surf.Pdata{1,2}.val = new_label';
Surf.Pdata{1,2}.name = 'new_map';
Surf.Face = Surf.Face-1;
vtkSurfWrite(output_surface_fname,Surf);

% inflate surface
infsurf = vtkSurfRead('/media/songyao/songyao/data/HCP_s900_data/display_surface/100206.inflated_MSMAll.32k_fs_LR.vtk');
infsurf.Face = infsurf.Face-1;
infsurf.Pdata=[];
infsurf.Pdata{1,1}.val = Surf.Pdata{1,1}.val;
infsurf.Pdata{1,1}.name = 'map';
infsurf.Pdata{1,2}.val = new_label;
infsurf.Pdata{1,2}.name = 'new_map';
vtkSurfWrite([surfdir,'groupwise_peaks_cluster_inf.vtk'],infsurf);

