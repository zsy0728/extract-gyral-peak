% all gyral peaks count map
clear
clc

output_dir = '/media/songyao/songyao/result/HCP_gyral_peak_MNINonLinear/groupwise_peaks/';
peakdir='/media/songyao/songyao/result/HCP_gyral_peak_MNINonLinear/individual_peaks_pits/';
surfdir = '/media/songyao/songyao/data/HCP_s900_data/T1_surface/';

tt='*.white_MSMAll.32k_fs_LR.vtk';
File = dir(fullfile(surfdir,tt));
FileNames = {File.name}';

sulc_all = zeros(length(FileNames),64984);
thick_all = zeros(length(FileNames),64984);
curv_all = zeros(length(FileNames),64984);
myelin_all = zeros(length(FileNames),64984);
countmap_peak = zeros(1,64984);

for sbj = 1:length(FileNames)
    disp(['sub = ',FileNames{sbj,1}(1:6)])
    surf = vtkSurfRead([surfdir,FileNames{sbj,1}(1:6),'.white_MSMAll.32k_fs_LR.vtk']);
    sulc = surf.Pdata{1,10}.val;
    curv = surf.Pdata{1,1}.val;
    myelin = surf.Pdata{1,7}.val;
    thick = surf.Pdata{1,11}.val;
    
    sulc_all(sbj,:) = sulc';
    thick_all(sbj,:) = thick';
    curv_all(sbj,:) = curv';
    myelin_all(sbj,:) = myelin';
    
    load([peakdir,FileNames{sbj,1}(1:6),'_peaks_ring=4.mat']);
    countmap_peak(local_maximum_vtxID_ring) = countmap_peak(local_maximum_vtxID_ring)+1;
end


sulc_all_mean = mean(sulc_all);
thick_all_mean = mean(thick_all);
curv_all_mean = mean(curv_all);
myelin_all_mean = mean(myelin_all);

infodir = '/media/songyao/a6eb3580-711b-4188-b899-9cc41af9e64d/songyao/result/HCP_gyral_peak_MNINonLinear/group_info/';
save([infodir,'S900_sulc_all.mat'],'sulc_all');
save([infodir,'S900_sulc_all_mean.mat'],'sulc_all_mean');
save([infodir,'S900_thick_all.mat'],'thick_all');
save([infodir,'S900_thick_all_mean.mat'],'thick_all_mean');
save([infodir,'S900_curv_all.mat'],'curv_all');
save([infodir,'S900_curv_all_mean.mat'],'curv_all_mean');
save([infodir,'S900_myelin_all.mat'],'myelin_all');
save([infodir,'S900_myelin_all_mean.mat'],'myelin_all_mean');

tempsurf = vtkSurfRead([surfdir,'100206.white_MSMAll.32k_fs_LR.vtk']);
tempsurf.Pdata=[];
tempsurf.Pdata{1,1}.val = countmap_peak;
tempsurf.Pdata{1,1}.name = 'peak_countmap';
tempsurf.Face = tempsurf.Face-1;
vtkSurfWrite([output_dir,'groupwise_peaks_count.vtk'],tempsurf);

tempsurf.Pdata=[];
tempsurf.Pdata{1,1}.val = sulc_all_mean;
tempsurf.Pdata{1,1}.name = 'sulc_all_mean';
tempsurf.Pdata{1,2}.val = thick_all_mean;
tempsurf.Pdata{1,2}.name = 'thick_all_mean';
tempsurf.Pdata{1,3}.val = curv_all_mean;
tempsurf.Pdata{1,3}.name = 'curv_all_mean';
tempsurf.Pdata{1,4}.val = myelin_all_mean;
tempsurf.Pdata{1,4}.name = 'myelin_all_mean';
vtkSurfWrite([infodir,'groupwise_features.vtk'],tempsurf);


%% smooth
clc
clear
addpath /media/songyao/a6eb3580-711b-4188-b899-9cc41af9e64d/songyao/matlab_path/Function_lin/
output_dir = '/media/songyao/a6eb3580-711b-4188-b899-9cc41af9e64d/songyao/result/HCP_gyral_peak_MNINonLinear/groupwise_peaks/';
surf_dir = '/media/songyao/a6eb3580-711b-4188-b899-9cc41af9e64d/songyao/data/HCP_s900_data/T1_surface/';

label_tag.POINT_DATA.SCALARS{1} = 'peak_countmap';
compsurf = ReadSurf([output_dir,'groupwise_peaks_count.vtk'],label_tag,0);
curvname = [surf_dir,'100206.white_MSMAll.32k_fs_LR_c.vtk.ply.vtk'];
cd(surf_dir)
vtk2ply_ascii([surf_dir,'100206.white_MSMAll.32k_fs_LR.vtk'],[surf_dir,'100206.white_MSMAll.32k_fs_LR_c.vtk.ply']);
system(['/media/songyao/a6eb3580-711b-4188-b899-9cc41af9e64d/songyao/matlab_path/ComputeCurvatures ',[surf_dir,'100206.white_MSMAll.32k_fs_LR_c.vtk.ply '],curvname]);
delete([surf_dir,'100206.white_MSMAll.32k_fs_LR_c.vtk.ply']);

surf_c = ReadSurf_2(curvname,{},1);
compsurf.Neighbors = surf_c.Neighbors;
label_tag.POINT_DATA.VECTORS{1} = 'Normals';
surf_normal = ReadSurf(curvname,label_tag,0);
normal = surf_normal.POINT_DATA.VECTORS{1};

ring_size = 1;
label = compsurf.POINT_DATA.SCALARS{1};
label_avg = zeros(1,size(surf_c.vertice,2));

for ii = 1:size(surf_c.vertice,2)
    den_ii = label(ii);
    center_ID = ii-1;
    [Neighbor_ID,Ring_ID] = Search_Neighbor_ID(compsurf,center_ID,ring_size,Inf);
    for iter=1:6
        sum1=0;
        sum2=0;
        for j=2:size(Neighbor_ID,2) %31
            num_j = Neighbor_ID(j)+1; %bianhao
            den_j = label(num_j);
            if normal(:,ii)'*normal(:,num_j)>0
                S = normal(:,ii)'*normal(:,num_j);
            else
                S=0;
            end
            dist = 1./(norm((compsurf.vertice(:,ii)-compsurf.vertice(:,num_j)),2));
            sum1 = sum1 + den_j*S*dist;
            sum2 = sum2 + S*dist;
        end
        if sum2==0
            ratio=0;
        else
            ratio = sum1/sum2;
        end
        
        den_ii=(den_ii+ratio)*0.5;
    end
    label_avg(ii) = den_ii;
end

outfname = [output_dir,'groupwise_peaks_count_smooth.vtk'];
label_tag3.POINT_DATA.SCALARS{1,1} = 'peaks_countmap_smooth';
compsurf.POINT_DATA.SCALARS{1,1} = label_avg;
WriteSurf(outfname,compsurf,label_tag3,'float',0);

% inflate surface
clc
clear
warning off

output_dir = '/media/songyao/songyao/result/HCP_gyral_peak_MNINonLinear/groupwise_peaks/';
surf1 = vtkSurfRead([output_dir,'groupwise_peaks_count.vtk']);
surf2 = vtkSurfRead([output_dir,'groupwise_peaks_count_smooth.vtk']);
inf_surf = vtkSurfRead(['/media/songyao/a6eb3580-711b-4188-b899-9cc41af9e64d/songyao/data/HCP_s900_data/display_surface/100206.inflated_MSMAll.32k_fs_LR.vtk']);
inf_surf.Pdata=[];
inf_surf.Pdata{1,1}.val = surf1.Pdata{1,1}.val;
inf_surf.Pdata{1,1}.name = 'peaks_countmap';
inf_surf.Pdata{1,2}.val = surf2.Pdata{1,1}.val;
inf_surf.Pdata{1,2}.name = 'peaks_countmap_smooth';
inf_surf.Face = inf_surf.Face-1;
vtkSurfWrite([output_dir,'groupwise_peaks_count_smooth_inf.vtk'],inf_surf);
