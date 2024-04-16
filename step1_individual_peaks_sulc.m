clear
clc

surf_dir = '/media/songyao/songyao/data/HCP_s900_data/T1_surface/';
output_dir = '/media/songyao/songyao/result/HCP_gyral_peak_MNINonLinear/individual_peaks_pits/';
tt='*.white_MSMAll.32k_fs_LR.vtk';
File = dir(fullfile(surf_dir,tt));
FileNames = {File.name}';
ring_size = 4;

for subid = 1:length(FileNames)
    sub = FileNames{subid}(1:6);
   
    surfname = [surf_dir,num2str(sub),'.white_MSMAll.32k_fs_LR.vtk'];
    if  ~exist(surfname,'file')
        continue;
    end
    disp(num2str(sub))
    surf = ReadSurf(surfname,[],0);
    surf2 = vtkSurfRead(surfname);
    sulc = surf2.Pdata{1,10}.val;
    thick = surf2.Pdata{1,11}.val;
    curv = surf2.Pdata{1,1}.val;
    curvname = [surf_dir,num2str(sub),'.white_MSMAll.32k_fs_LR_c.vtk'];
    surf_c = ReadSurf_2(curvname,{},1);
    surf.Neighbors = surf_c.Neighbors;
    
    %single point gyral peaks(sulc)
%     local_maximum_vtxID_sulc=[];
%     indx=1;
%     for i = 1:size(surf.vertice,2)
%         if curv(i)<0||thick(i)==0
%             continue;
%         end
%         [Neighbor_ID,Ring_ID] = Search_Neighbor_ID(surf,i-1,ring_size,Inf);
%         tmp_label = sulc(Neighbor_ID+1);
%         if max(tmp_label) == tmp_label(1)
%             local_maximum_vtxID_sulc(indx) = i;
%             indx = indx + 1;
%         else
%         end
%     end
%     vtx={};
%     outfname = [output_dir,num2str(sub),'_single_peak_ring=',num2str(ring_size),'.vtk'];
%     vtx.vertice = surf.vertice(:,local_maximum_vtxID_sulc);
%     VertexWrite(outfname,vtx,{},'float');
%     save([output_dir,num2str(sub),'_single_peak_ring=',num2str(ring_size),'.mat'],'local_maximum_vtxID_sulc');
    
    %ring gyral peaks
    indx=1;
    local_maximum_vtxID_ring=[];
    ring_size2 = 1;
    for i = 1:size(surf.vertice,2)
        if curv(i)<0||thick(i)==0
            continue;
        end
        [Neighbor_ID,Ring_ID] = Search_Neighbor_ID(surf_c,i-1,ring_size,Inf);
        tmp_label = sulc(Neighbor_ID+1);
        if max(tmp_label) == tmp_label(1)
            [Max_Neighbor_ID,Max_Ring_ID] = Search_Neighbor_ID(surf_c,i-1,ring_size2,Inf);
            local_maximum_vtxID_ring = [local_maximum_vtxID_ring,(Max_Neighbor_ID+1)];
            indx = indx + 1;
        else
        end
    end
    save([output_dir,num2str(sub),'_peaks_ring=',num2str(ring_size),'.mat'],'local_maximum_vtxID_ring');
    outfname = [output_dir,num2str(sub),'_peaks_ring=',num2str(ring_size),'.vtk'];
    vtx.vertice = surf.vertice(:,local_maximum_vtxID_ring);
    VertexWrite(outfname,vtx,{},'float');
end

% inflate surface
% clc
% clear
% warning off
% 
% output_dir = '/media/songyao/songyao/result/HCP_gyal_peak/individual_peaks_sulc/';
% load([output_dir,'100206_peaks_ring=4.mat'])
% inf_surf = vtkSurfRead('/media/songyao/songyao/data/HCP_s900_data/display_surface/100206.inflated_MSMAll.32k_fs_LR.vtk');
% outfname = [output_dir,'100206_peaks_ring=4_inf.vtk'];
% vtx.vertice = inf_surf.Vtx(:,local_maximum_vtxID_ring);
% VertexWrite(outfname,vtx,{},'float');
% 
% load([output_dir,'100206_single_peak_ring=4.mat'])
% outfname = [output_dir,'100206_single_peak_ring=4_inf.vtk'];
% vtx.vertice = inf_surf.Vtx(:,local_maximum_vtxID_sulc);
% VertexWrite(outfname,vtx,{},'float');


