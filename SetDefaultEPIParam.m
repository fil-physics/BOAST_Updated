function epi_param = SetDefaultEPIParam

% =========================================================================
% This function sets the fixed parameters for the EPI protocol. 
% These values can be modified as needed to suit specific requirements 
% or preferences.
% =========================================================================
% main_orientation                : Slice oriantation
%                                   'TRA' : transverse 
%                                   'CRO' : coronal
%                                   'SAG' : sagittal
% fov                             : Field of view (in mm)
% base_res                        : Basic resolution (Matrix size)
% pe_ov                           : Oversampling Ratio in Phase Encoding 
%                                   Direction in %
% PF                              : Partial Fourier Coefficient
% AF                              : In-plane Acceleration Factor
% slicethickness                  : Full width at half-maximum (FWHM) 
%                                   of the slice profile (in mm)
% echo_spacing                    : Echo spacing (in ms)
% echotime                        : Effective (central) echo time (in ms)
% vox                             : Voxel size (in mm) 
%                                   (1x3) array (read, phase, slice) direction
% =========================================================================

% Updated 23/09/2024
% by Shokoufeh Golshani

epi_param.main_orientation = 'TRA';
epi_param.fov              = 192;
epi_param.base_res         = 64;
epi_param.pe_ov            = 13;
epi_param.PF               = 1;
epi_param.AF               = 1;
% Note here 2 mm is used as the FWHM
epi_param.slicethickness   = 2; % Siemens pulse approximates Gaussian with 2 mm FWHM
epi_param.echo_spacing     = 0.5; 
epi_param.echotime= 30; 
epi_param.vx               = [3 3 3];



end
