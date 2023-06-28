%% Load the data
subj1 = load('Subj1_Betas_288sounds.mat');
subj2 = load('Subj2_Betas_288sounds.mat');
subj3 = load('Subj3_Betas_288sounds.mat');
subj4 = load('Subj4_Betas_288sounds.mat');
subj5 = load('Subj5_Betas_288sounds.mat');

%% Concatenate the fMRI data for all subjects and regions
n_subjects = 5;
n_regions = 4;

% Initialize the concatenated matrices
concatenatedMatrices = cell(n_regions, 1);

for subject = 1:n_subjects
    % Get the subject structure for the current subject
    subjectData = eval(sprintf('subj%d', subject));
    
    for region = 1:n_regions
        % Get the voxel data for the current brain region
        regionData = subjectData.Data(region).SoundResponse;
        
        % Concatenate the voxel data for the current brain region
        if isempty(concatenatedMatrices{region})
            concatenatedMatrices{region} = regionData;
        else
            concatenatedMatrices{region} = horzcat(concatenatedMatrices{region}, regionData);

        end
    end
end

%% Check that the number of voxels for all the subjects is correct
n_vox_A1_allsubj = length(subj1.Data(1).SoundResponse) + length(subj2.Data(1).SoundResponse) + length(subj3.Data(1).SoundResponse) + length(subj4.Data(1).SoundResponse) + length(subj5.Data(1).SoundResponse);
n_vox_R_allsubj = length(subj1.Data(2).SoundResponse) + length(subj2.Data(2).SoundResponse) + length(subj3.Data(2).SoundResponse) + length(subj4.Data(2).SoundResponse) + length(subj5.Data(2).SoundResponse);
n_vox_S_allsubj = length(subj1.Data(3).SoundResponse) + length(subj2.Data(3).SoundResponse) + length(subj3.Data(3).SoundResponse) + length(subj4.Data(3).SoundResponse) + length(subj5.Data(3).SoundResponse);
n_vox_F_allsubj = length(subj1.Data(4).SoundResponse) + length(subj2.Data(4).SoundResponse) + length(subj3.Data(4).SoundResponse) + length(subj4.Data(4).SoundResponse) + length(subj5.Data(4).SoundResponse);

n_vox_A1_allsubj == length(concatenatedMatrices{1})
n_vox_R_allsubj == length(concatenatedMatrices{2})
n_vox_S_allsubj == length(concatenatedMatrices{3})
n_vox_F_allsubj == length(concatenatedMatrices{4})

% All are correct, so I will store the values in a matrix to export it as a
% csv file for having a reference
voxels_regions_idx = [n_vox_A1_allsubj, n_vox_R_allsubj, n_vox_S_allsubj, n_vox_F_allsubj];

voxels_subj1_idx = [length(subj1.Data(1).SoundResponse), length(subj1.Data(2).SoundResponse), length(subj1.Data(3).SoundResponse), length(subj1.Data(4).SoundResponse)];

voxels_subj2_idx = [length(subj2.Data(1).SoundResponse), length(subj2.Data(2).SoundResponse), length(subj2.Data(3).SoundResponse), length(subj2.Data(4).SoundResponse)];

voxels_subj3_idx = [length(subj3.Data(1).SoundResponse), length(subj3.Data(2).SoundResponse), length(subj3.Data(3).SoundResponse), length(subj3.Data(4).SoundResponse)];

voxels_subj4_idx = [length(subj4.Data(1).SoundResponse), length(subj4.Data(2).SoundResponse), length(subj4.Data(3).SoundResponse), length(subj4.Data(4).SoundResponse)];

voxels_subj5_idx = [length(subj5.Data(1).SoundResponse), length(subj5.Data(2).SoundResponse), length(subj5.Data(3).SoundResponse), length(subj5.Data(4).SoundResponse)];

% Create a cell array with the variable names
regions = {'A1', 'R', 'Slow', 'Fast'};

% Create a table with the variable names as a column and the array values
voxels_regions_idx_table = table(regions', voxels_regions_idx', 'VariableNames', {'Region', 'Allsubjects'});

voxels_regions_idx_table.Subject1 = voxels_subj1_idx';
voxels_regions_idx_table.Subject2 = voxels_subj2_idx';
voxels_regions_idx_table.Subject3 = voxels_subj3_idx';
voxels_regions_idx_table.Subject4 = voxels_subj4_idx';
voxels_regions_idx_table.Subject5 = voxels_subj5_idx';
%% Combine all regions for al subjects in a single matrix
fmri_allregions_and_subj = horzcat(concatenatedMatrices{:});

%% Export concatenated fmri data and number of voxels per region
filename_fmri_allregions = 'fmri_allsubj_and_regions.csv';
writematrix(fmri_allregions_and_subj, filename_fmri_allregions);

filename_voxel_regions_idx = 'fmri_voxel_regions_idx.csv';
writetable(voxels_regions_idx_table, filename_voxel_regions_idx);

filename_fmri_A1 = 'fmri_allsubj_A1.csv';
writematrix(concatenatedMatrices{1}, filename_fmri_A1);

filename_fmri_R = 'fmri_allsubj_R.csv';
writematrix(concatenatedMatrices{2}, filename_fmri_R);

filename_fmri_S = 'fmri_allsubj_S.csv';
writematrix(concatenatedMatrices{3}, filename_fmri_S);

filename_fmri_F = 'fmri_allsubj_F.csv';
writematrix(concatenatedMatrices{4}, filename_fmri_F);


